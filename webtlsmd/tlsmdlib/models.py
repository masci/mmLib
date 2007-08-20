from django.db import models
from django.contrib.auth.models import User
from django.conf import settings

import os



class TlsmdJobException(Exception):
    pass


class TlsmdJob(models.Model):
    STATE_CHOICES = (
        ('submit1', 'Submit Step 1'),
        ('submit2', 'Submit Step 2'),
        ('queued', 'Job Queued'),
        ('completed', 'Job Completed'),
        ('running', 'Job Running'),
        ('failed', 'Job Failed'),
        ('removed', 'Job Removed'),
        )
    state = models.CharField(maxlength=32, choices=STATE_CHOICES, default='submit1')

    plot_format = models.CharField(maxlength=32, choices=(('PNG','PNG'),('SVG','SVG')), default='PNG')
    tls_model = models.CharField(maxlength=32, choices=(('ISOT', 'Isotropic'), ('ANISO', 'Anisotropic')), default='ISOT')
    include_atoms = models.CharField(maxlength=32, choices=(('ALL', 'All Atoms'), ('MAINCHAIN', 'Main Chain Atoms')), default='ALL')
    weight_function = models.CharField(maxlength=32, choices=(('NONE', 'None'), ('IUISO', 'Inverse Uiso Magnitude')), default='NONE')

    user = models.ForeignKey(User, null=True, blank=True)
    private_job = models.BooleanField(default=False)
    full_name = models.CharField(maxlength=255, null=True, blank=True)
    email = models.CharField(maxlength=255, null=True, blank=True)
    ip_address = models.CharField(maxlength=32, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)

    structure_id = models.CharField(maxlength=32, null=True, blank=True, db_index=True)
    pdb_filename = models.CharField(maxlength=255, default='struct.pdb')
    log_filename = models.CharField(maxlength=255, default='log.txt')
    analysis_dirname = models.CharField(maxlength=255, default='ANALYSIS')

    created = models.DateTimeField(auto_now_add=True)
    last_modified = models.DateTimeField(auto_now_add=True, auto_now=True)
    run_time_begin = models.DateTimeField(null=True, blank=True)
    run_time_end = models.DateTimeField(null=True, blank=True)

    class Meta:
        db_table = 'tlsmd_job'

    def get_job_dirname(self):
        assert self.id is not None
        return 'TLSMD%s' % self.id
    job_dirname = property(get_job_dirname)

    def get_job_dir(self):
        return os.path.join(settings.TLSMD_WORK_DIR, self.get_job_dirname())
    job_dir = property(get_job_dir)

    def get_structure_path(self):
        return os.path.join(self.get_job_dir(), self.pdb_filename)

    def get_job_url(self):
        return '/tlsmd/job/%d/' % self.id
    job_url = property(get_job_url)

    def get_job_output_baseurl(self):
        return '%s/%s' % (settings.TLSMD_WORK_URL, self.get_job_dirname())
    job_output_baseurl = property(get_job_output_baseurl)

    def get_log_path(self):
        return os.path.join(self.get_job_dir(), self.log_filename)
    log_path = property(get_log_path)

    def get_log_url(self):
        if not os.path.isfile(self.get_log_path()):
            return None
        return '%s/%s' % (self.get_job_output_baseurl(), self.log_filename)
    log_url = property(get_log_url)

    def get_analysis_dir(self):
        return os.path.join(self.get_job_dir(), self.analysis_dirname)
    analysis_dir = property(get_analysis_dir)

    def get_analysis_url(self):
        if not os.path.isfile(os.path.join(self.get_analysis_dir, 'index.html')):
            return None
        return '%s/%s/index.html' % (self.get_job_output_baseurl(), self.analysis_dirname)
    analysis_url = property(get_analysis_url)

    def set_structure_file(self, struct_file_data):
        """
        Creates job directory, saves structure file to the job directory, sets
        this job up to run.
        """
        if not os.path.isdir(settings.TLSMD_WORK_DIR):
            raise TlsmdJobException('settings.TLSMD_WORK_DIR=%s directory does not exist' % settings.TLSMD_WORK_DIR)

        if os.path.exists(self.get_job_dir()):
            self.remove_job_dir()
            
        try:
            os.mkdir(self.get_job_dir())
        except OSError:
            raise TlsmdJobException('cannot make directory=%s' % self.get_job_dir())

        try:
            open(self.get_structure_path(),'w').write(struct_file_data)
        except IOError:
            raise TlsmdJobException('cannot write file=%s' % self.get_structure_path())

        ## load the structure with mmLib, extract data from the structure
        from mmLib import FileIO
        try:
            struct = FileIO.LoadStructure(fil = self.get_structure_path())
        except:
            raise TlsmdJobException('The Python Macromolecular Library was unable to load your structure file.')

        ## set self.structure_id if necessary
        if not self.structure_id and struct.structure_id:
            self.structure_id = struct.structure_id
        if not self.structure_id:
            self.structure_id = 'XXXX'

        ## Select Chains for Analysis
        num_atoms = 0
        num_aniso_atoms = 0
        largest_chain_seen = 0

        ## catelog all the chains of the structure
        chains = []
        for chain in struct.iter_chains():
            job_chain = TlsmdJobChain(tlsmd_job=self, chain_id=chain.chain_id)
            job_chain.set_chain(chain)
            num_atoms += job_chain.num_atoms
            num_aniso_atoms += job_chain.num_aniso_atoms
            largest_chain_seen = max(job_chain.length, largest_chain_seen)

        if num_atoms < 1:
            self.remove_job()
            raise TlsmdJobException('Your submitted structure contained no atoms')

        if largest_chain_seen > settings.MAX_CHAIN_LENGTH:
            self.remove_job()
            raise TlsmdJobException('Your submitted structure contained a chain exceeding the %d residue limit' % settings.MAX_CHAIN_LENGTH)

        try:
            aniso_ratio = float(num_aniso_atoms) / float(num_atoms)
        except ZeroDivisionError:
            raise TlsmdJobException('Your submitted structure contained no atoms')

        if aniso_ratio > 0.90:
            self.tls_model = 'ANISO'

        self.save()

    def remove_job_dir(self):
        job_dir = self.get_job_dir()
        if job_dir and job_dir.startswith(settings.TLSMD_WORK_DIR) and os.path.isdir(job_dir):
            for root, dirs, files in os.walk(job_dir, topdown = False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
            os.rmdir(job_dir)

    def remove_job(self):
        """
        Removes all the files for a job and marks the job state as removed.
        """
        self.remove_job_dir()
        self.state = 'removed'
        self.save()

    def possible_chains(self):
        """
        Returns a list of all selected chains.
        """
        return list(self.chains.filter(tlsmd_analysis_possible=True))

    def selected_chains(self):
        """
        Returns a list of all selected chains.
        """
        return list(self.chains.filter(selected=True, tlsmd_analysis_possible=True))


            

class TlsmdJobChain(models.Model):
    tlsmd_job = models.ForeignKey('TlsmdJob', related_name='chains', db_index=True)
    chain_id = models.CharField(maxlength=32)
    name = models.CharField(maxlength=32)
    length = models.IntegerField()
    num_atoms = models.IntegerField(default=0)
    num_aniso_atoms = models.IntegerField(default=0)
    preview = models.CharField(maxlength=255)
    description = models.CharField(maxlength=255)
    tlsmd_analysis_possible = models.BooleanField()
    selected = models.BooleanField()
    
    class Meta:
        db_table = 'tlsmd_job_chain'
        unique_together = [('tlsmd_job', 'chain_id')]

    def set_chain(self, chain):
        assert self.chain_id == chain.chain_id

        self.name = 'CHAIN%s' % (self.chain_id)

        ## number of fragments
        naa = chain.count_amino_acids()
        nna = chain.count_nucleic_acids()
        num_frags = 0
        if naa > 0:
            num_frags = naa
        elif nna > 0:
            num_frags = nna
        self.length = num_frags

        ## less than 10 fragments? don't select
        if num_frags < 10:
            self.tlsmd_analysis_possible = False
            self.selected = False
        else:
            self.tlsmd_analysis_possible = True
            self.selected = True 

        ## create chain description 
        if naa > 0:
            self.description = 'Chain %s (%d Amino Acid Residues)' % (self.chain_id, num_frags)
        elif nna > 0:
            self.description = 'Chain %s (%d Nucleic Acid Residues)' % (self.chain_id, num_frags)
        else:
            self.description = 'Chain %s' % self.chain_id

        num_atoms = 0
        num_aniso_atoms = 0
        for atm in chain.iter_all_atoms():
            num_atoms += 1
            if atm.U is not None:
                num_aniso_atoms += 1
        self.num_atoms = num_atoms
        self.num_aniso_atoms = num_aniso_atoms

        ## chain sequence preview
        listx = []
        i = 0
        for frag in chain.iter_fragments():
            i += 1
            if i > 5:
                break
            listx.append(frag.res_name)
        self.preview = ' '.join(listx)

        self.save()
