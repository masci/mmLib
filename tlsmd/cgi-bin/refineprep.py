## TLS Minimized Domains (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import time
import socket
import string
import glob

import xmlrpclib
import cgitb; cgitb.enable()
import cgi

from mmLib.Structure      import *
from mmLib.FileLoader     import *
from mmLib.Extensions.TLS import *

## GLOBALS
from cgiconfig import *
webtlsmdd = xmlrpclib.ServerProxy(WEBTLSMDD, allow_none=True)

CAPTION = """\
Download both this modified PDB file of your structure, and the corresponding
TLSIN file. Feed these to REFMAC5 as a start point for  multi-TLS group refinement.
See the TLSMD documentation for detailed instructions.
"""

class UInvalid(Exception):
    """Exception raised if a anisotropic ADP U is not positive definite.
    """
    def __init__(self, atom, U):
        Exception.__init__(self)
        self.text = "%s eigenvalues(%s)" % (str(atom), str(eigenvalues(U2B*U))) 

    def __str__(self):
        return self.text


def refmac5_prep(xyzin, tlsin_list, xyzout, tlsout):
    """Use TLS model + Uiso for each atom.  Output xyzout with the
    residual Uiso only.
    """
    os.umask(022)
    
    ## load structure
    struct = LoadStructure(fil = xyzin)

    ## load and construct TLS groups
    tls_group_list = []

    tls_file = TLSFile()
    tls_file.set_file_format(TLSFileFormatTLSOUT())

    tls_file_format = TLSFileFormatTLSOUT()

    for tlsin in tlsin_list:
        tls_desc_list = tls_file_format.load(open(tlsin, "r"))

        for tls_desc in tls_desc_list:
            tls_file.tls_desc_list.append(tls_desc)
            tls_group = tls_desc.construct_tls_group_with_atoms(struct)
            tls_group.tls_desc = tls_desc
            tls_group_list.append(tls_group)

    ## set the extra Uiso for each atom
    for tls_group in tls_group_list:

        ## minimal/maximal amount of Uiso which has to be added
        ## to the group's atoms to to make Uiso == Uiso_tls
        min_Uiso = 0.0
        max_Uiso = 0.0

        n         = 0
        sum_diff2 = 0.0

        for atm, Utls in tls_group.iter_atm_Utls():
            for aatm in atm.iter_alt_loc():
                tls_tf = trace(Utls)/3.0
                ref_tf = trace(aatm.get_U())/3.0

                n += 1
                sum_diff2 += (tls_tf - ref_tf)**2
                
                if ref_tf>tls_tf:
                    max_Uiso = max(ref_tf - tls_tf, max_Uiso)
                else:
                    min_Uiso = max(tls_tf - ref_tf, min_Uiso)

        msd = sum_diff2 / n
        rmsd = math.sqrt(msd)


        ## report the percentage of atoms with Uiso within the RMSD
        ntotal = 0
        nrmsd  = 0
        
        for atm, Utls in tls_group.iter_atm_Utls():
            for aatm in atm.iter_alt_loc():
                tls_tf = trace(Utls)/3.0
                ref_tf = trace(aatm.get_U())/3.0

                ntotal += 1
                deviation = math.sqrt((tls_tf - ref_tf)**2)
                
                if deviation<=rmsd:
                    nrmsd += 1

        ## reduce the TLS group T tensor by min_Uiso so that
        ## a PDB file can be written out where all atoms
        ## Uiso == Uiso_tls

        ## we must rotate the T tensor to its primary axes before
        ## subtracting min_Uiso magnitude from it
        (T_eval, TR) = eigenvectors(tls_group.T)
        T = matrixmultiply(TR, matrixmultiply(tls_group.T, transpose(TR)))

        assert allclose(T[0,1], 0.0)
        assert allclose(T[0,2], 0.0)
        assert allclose(T[1,2], 0.0)

        T[0,0] = T[0,0] - min_Uiso
        T[1,1] = T[1,1] - min_Uiso
        T[2,2] = T[2,2] - min_Uiso

        ## now take half of the smallest principal component of T and
        ## move it into the individual atomic temperature factors

        min_T    = min(T[0,0], min(T[1,1], T[2,2]))
        sub_T    = min_T * 0.80
        add_Uiso = min_T - sub_T
        
        T[0,0] = T[0,0] - sub_T
        T[1,1] = T[1,1] - sub_T
        T[2,2] = T[2,2] - sub_T
        
        ## rotate T back to original orientation
        tls_group.T = matrixmultiply(transpose(TR), matrixmultiply(T, TR))

        ## reset the TLS tensor values in the TLSDesc object so they can be saved
        tls_group.tls_desc.set_tls_group(tls_group)
        
        ## set atm.temp_factor
        for atm, Utls in tls_group.iter_atm_Utls():
            for aatm in atm.iter_alt_loc():
                tls_tf = trace(Utls)/3.0
                ref_tf = trace(aatm.get_U())/3.0
                
                if ref_tf>tls_tf:
                    aatm.temp_factor = ((add_Uiso) + ref_tf - tls_tf)*U2B
                    aatm.U = None
                else:
                    aatm.temp_factor = (add_Uiso) * U2B
                    aatm.U = None

    SaveStructure(fil=xyzout, struct=struct)
    tls_file.save(open(tlsout, "w"))
    


class Page(object):
    def __init__(self, form):
        self.form = form

    def html_title(self, title):
        return '<center><h1>%s</h1></center>' % (title)

    def html_head_nocgi(self, title):
        x  = ''
        x += '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" '
        x += '"http://www.w3.org/TR/html4/loose.dtd">\n\n'
        x += '<html>'
        x += '<head>'
        x += '  <title>%s</title>' % (title)
        x += '  <style type="text/css" media=screen>'
        x += '  <!-- '
        x += '  BODY { background-color: white;'
        x += '         margin-left: 5%; margin-right: 5%;'
        x += '         border-left: 5%; border-right: 5%;'
        x += '         margin-top: 2%; border-top: 2%;}'
        x += '  -->'
        x += '  </style>'
        x += '</head>'
        x += '<body>'
        return x

    def html_head(self, title):
        x = ''
        x += 'Content-Type: text/html\n\n'
        x += self.html_head_nocgi(title)
        return x

    def html_foot(self):
        x = ''
        x += '<center>'
        x += '<p><small><b>Version %s</b> Last Modified %s' % (VERSION, LAST_MODIFIED_DATE)
        x += ' by %s ' % (LAST_MODIFIED_BY)
        x += '<i>%s</i></small></p>' % (LAST_MODIFIED_BY_EMAIL)
        x += '</center>'
        x += '</body></html>'
        return x


class ErrorPage(Page):
    def __init__(self, form, text):
        Page.__init__(self, form)
        self.text = text
    
    def html_page(self):
        title = 'TLSMD: Error'
        
        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)
        x += '<br>'
        x += '<center><h3>A Error Occured</h3></center>'

        if self.text!=None:
            x += self.text
            
        x += self.html_foot()
        return x


class RefinePrepError(Exception):
    def __init__(self, text):
        Exception.__init__(self)
        self.text = text


class RefinePrepPage(Page):
    def html_page(self):
        job_id = check_job_id(self.form)
        
        title = 'Input Files for Refmac5 TLS Refinement'
        
        x  = ''
        x += self.html_head(title)
        x += self.html_title(title)

        x += '<center>'
        x += '<h3>'
        x += 'Step 2: Download the generated XYZIN(PDBIN) and TLSIN files below'
        x += '</h3>'
        x += '</center>'
        
        ## get the analysis directory, and make sure it exists
        analysis_dir = webtlsmdd.job_data_get(job_id, "analysis_dir")
        if not os.path.isdir(analysis_dir):
            raise RefinePrepError("Analysis Directory Not Found")
        os.chdir(analysis_dir)
        
        ## extract ntls selections from CGI form
        chain_ntls = []
        for key in self.form.iterkeys():
            if key.startswith("NTLS_CHAIN"):
                chain_id = key[-1]
                try:
                    ntls = int(self.form[key].value)
                except ValueError:
                    continue
                chain_ntls.append((chain_id, ntls))

        chain_ntls.sort()

        ## make sure there were selections
        if len(chain_ntls)==0:
            raise RefinePrepError("Form Processing Error: No Chains Selected")

        ## now create filenames
        struct_id = webtlsmdd.job_data_get(job_id, "structure_id")
        assert struct_id!=False
        assert struct_id

        ## input structure
        pdbin  = "%s.pdb" % (struct_id)
        if not os.path.isfile(pdbin):
            pdbin = None
            for pdbx in glob.glob("*.pdb"):
                if len(pdbx)==8:
                    struct_id = pdbx[:4]
                    pdbin = pdbx
                    break
            if pdbin==None:	
                raise RefinePrepError("Input PDB File %s Not Found" % (pdbin))

        ## the per-chain TLSOUT files from TLSMD must be merged
        tlsins = []
        for chain_id, ntls in chain_ntls:
            tlsin = "%s_CHAIN%s_NTLS%d.tlsout" % (struct_id, chain_id, ntls)
            if not os.path.isfile(tlsin):
                raise RefinePrepError("Input TLSIN File %s Not Found" % (tlsin))
            tlsins.append(tlsin)

        ## form unique pdbout/tlsout filenames
        listx = [struct_id]
        for chain_id, ntls in chain_ntls:
            listx.append("CHAIN%s" % (chain_id))
            listx.append("NTLS%d" % (ntls))
        outbase = string.join(listx, "_")
        pdbout = "%s.pdb" % (outbase)

        ## the tlsout from this program is going to be the tlsin
        ## for refinement, so it's important for the filename to have
        ## the tlsin extension so the user is not confused
        tlsout = "%s.tlsin" % (outbase)

        ## make urls for linking
        analysis_base_url = webtlsmdd.job_data_get(job_id, "analysis_base_url")
        pdbout_url = "%s/%s" % (analysis_base_url, pdbout)
        tlsout_url = "%s/%s" % (analysis_base_url, tlsout)

        ## create the files
        refmac5_prep(pdbin, tlsins, pdbout, tlsout)

        x += '<p>%s</p>' % (CAPTION)

        ## success -- make download links
        x += '<center>'
        x += '<table border="0" style="background-color:#eeeeee">'
        x += '<tr>'
        x += '<td align="right"><b>PDBIN File</b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (pdbout_url, pdbout)
        x += '</tr><tr>'
        x += '<td align="right"><b>TLSIN File</b></td>'
        x += '<td><a href="%s" type="text/plain">%s</a></td>' % (tlsout_url, tlsout)
        x += '</table>'

        x += '<br>'

        x += '<center>'
        x += '<h3>'
        x += 'Step 3: Read this '
        x += '<a href="/~tlsmd/documentation.html#refmac5">How-To</a>'
        x += '</h3>'
        x += '</center>'

        x += self.html_foot()
        return x


def check_job_id(form):
    """Retrieves and confirms the job_id from a incomming form.  Returns
    None on error, or the job_id on success.
    """
    if form.has_key("job_id"):
        job_id = form["job_id"].value
        if len(job_id)<20:
            if job_id.startswith("TLSMD"):
                if webtlsmdd.job_exists(job_id):
                    return job_id
    return None


def main():
    form = cgi.FieldStorage()

    page = None
    job_id = check_job_id(form)
    if job_id==None:
        page = ErrorPage(form, "<center>The Job ID seems to be expired or invalid.</center>")
    else:
        page = RefinePrepPage(form)


    try:
        print page.html_page()

    except RefinePrepError, err:
        text = '<center><p>%s</p></center>' % (err.text)
        page = ErrorPage(form, text)
        print page.html_page()

    except xmlrpclib.Fault, err:
        page = ErrorPage(form, "xmlrpclib.Fault: " +str(err))
        print page.html_page()

    except socket.error, err:
        page = ErrorPage(form, "socket.error: " + str(err))
        print page.html_page()


if __name__=="__main__":
    main()
    sys.exit(0)
