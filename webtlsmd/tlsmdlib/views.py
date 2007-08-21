import os
import sys
import time
import re
import string
import random

from django.conf import settings
from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext
from django.shortcuts import render_to_response
from django.contrib.auth.decorators import login_required


##
## Views
##

def default_redirect(request):
    return HttpResponseRedirect('/tlsmd/')



def index(request):
    return render_to_response('webtlsmd/index.html', RequestContext(request))



def submit(request):
    """
    First-stage job submission.
    The user uploads a PDB/mmCIF file, or enters a PDB ID for the the system to auto-download
    from the RCSB.
    """
    if isinstance(request.session.get('submit_job_id'), (int, long)):
        return submit_job(request, request.session['submit_job_id'])

    if request.method == 'POST':

        ## uploaded file
        if request.FILES.has_key('pdbfile'):

            from tlsmdlib.models import TlsmdJob, TlsmdJobException
            tlsmd_job = TlsmdJob.objects.create()
            if not request.user.is_anonymous():
                tlsmd_job.user = request.user
            tlsmd_job.ip_address = request.META.get('REMOTE_ADDR', 'Unknown')
            tlsmd_job.set_structure_file(request.FILES['pdbfile']['content'])

            request.session['submit_job_id'] = tlsmd_job.id

        return HttpResponseRedirect(request.path)

    rc = RequestContext(request)
    return render_to_response('webtlsmd/submit1.html', rc)



def submit_job(request, job_id):
    """
    Second-stage job submission.
    The user selects various TLSMD options before final submission.
    """
    from tlsmdlib.models import TlsmdJob
    try:
        tlsmd_job = TlsmdJob.objects.get(id=job_id)
    except TlsmdJob.DoesNotExist:
        del request.session['submit_job_id']
        return HttpResponseRedirect(request.path)

    ## check the user is the same
    if tlsmd_job.user is not None:
        assert request.user.id == tlsmd_job.user_id

    if request.method == 'POST':
        form = request.POST

        ## check for job submission cancel 
        if form.get('cancel'):
            del request.session['submit_job_id']
            tlsmd_job.remove_job()
            return HttpResponseRedirect(request.path)

        form_valid = True

        if form.has_key("private_job"):
            tlsmd_job.private_job = True

        if form.has_key("user_name"):
            tlsmd_job.full_name = form["user_name"].strip()
        else:
            form_valid = False

        if form.has_key("email"):
            email_address = form["email"].strip()
            if len(email_address) > 5:
                tlsmd_job.email = email_address
            else:
                form_valid = False

        if form.has_key("structure_id"):
            structure_id = form["structure_id"].strip()
            if vet_data(structure_id, 4):
                tlsmd_job.structure_id = structure_id

        for job_chain in tlsmd_job.chains.all():
            if form.has_key(job_chain.name):
                job_chain.selected = True
            else:
                job_chain.selected = False
            job_chain.save()

        if form.has_key("tls_model"):
            tls_model = form["tls_model"].strip()
            if tls_model in ["ISOT", "ANISO"]:
                tlsmd_job.tls_model = tls_model

        if form.has_key("weight"):
            weight = form["weight"].strip()
            if weight in ["NONE", "IUISO"]:
                tlsmd_job.weight_function = weight

        if form.has_key("include_atoms"):
            include_atoms = form["include_atoms"].strip()
            if include_atoms in ["ALL", "MAINCHAIN"]:
                tlsmd_job.include_atoms = include_atoms

        if form.has_key("plot_format"):
            plot_format = form["plot_format"].strip()
            if plot_format in ["PNG", "SVG"]:
                tlsmd_job.plot_format = plot_format

        ## keep self re-directing until the form is valid
        if not form_valid:
            tlsmd_job.save()
            return HttpResponseRedirect(request.path)

        ## if everything with the form is okay, then change
        ## the job state to queued
        del request.session['submit_job_id']
        
        tlsmd_job.state = 'queued'
        tlsmd_job.save()
        return HttpResponseRedirect(tlsmd_job.job_url)


    tc = {'job': tlsmd_job}
    rc = RequestContext(request, tc)
    return render_to_response('webtlsmd/submit2.html', rc)


    
    
SELECTED_JOBS_RE = re.compile('^select_(\d+)$')
def get_selected_jobs(kwvars):
    selected_jobs = []
    for key in kwvars.keys():
        mobj = SELECTED_JOBS_RE.match(key)
        if mobj:
            selected_jobs.append(int(mobj.group(1)))
    return selected_jobs
    

def jobs(request):
    from tlsmdlib.models import TlsmdJob

    if request.method == 'POST':
        selected_jobs = get_selected_jobs(request.POST)
        if len(selected_jobs) > 0:
            for job_id in selected_jobs:
                try:
                    tlsmd_job = TlsmdJob.objects.get(id=job_id)
                except TlsmdJob.DoesNotExist:
                    continue
                else:
                    tlsmd_job.remove_job()
        
        return HttpResponseRedirect(request.path)


    running_job_list = list(TlsmdJob.objects.filter(state='running').order_by('id'))
    completed_job_list = list(TlsmdJob.objects.filter(state__in=['running','error']).order_by('id'))
    queued_job_list = list(TlsmdJob.objects.filter(state='queued').order_by('id'))
    limbo_job_list = list(TlsmdJob.objects.all().exclude(state__in=['running','error','queued']).order_by('id'))

    job_lists = [
        {'title': 'Running Jobs', 'job_list': running_job_list},
        {'title': 'Queued Jobs',  'job_list': queued_job_list},
        {'title': 'Completed Jobs', 'job_list': completed_job_list},
        {'title': 'Broken Jobs', 'job_list': limbo_job_list} ]
        
    tc = {
        'job_lists': job_lists,
        'administrator': True }

    rc = RequestContext(request, tc)
    return render_to_response('webtlsmd/jobs.html', rc)


def job(request, job_id):
    job_id = int(job_id)
    
    from tlsmdlib.models import TlsmdJob
    try:
        tlsmd_job = TlsmdJob.objects.get(id=job_id)
    except TlsmdJob.DoesNotExist:
        return HttpResponse('Sorry, the job %d does not exist' % job_id)

    if tlsmd_job.user and tlsmd_job.private_job:
        if request.user.id != tlsmd_job.user_id:
            return HttpResponse('Sorry, you are not authorized to view this structure.')

    rc = RequestContext(request, {'job':tlsmd_job, 'debug_vars':False})
    return render_to_response('webtlsmd/job.html', rc)


def examples(request):
    return render_to_response('webtlsmd/examples.html', RequestContext(request))


def documentation(request):
    return render_to_response('webtlsmd/documentation.html', RequestContext(request))


@login_required
def my_account(request):
    return render_to_response('webtlsmd/my_account.html', RequestContext(request))


def vet_data(data, max_len):
    if isinstance(data, unicode):
        return False
    if len(data) > max_len:
        return False
    if not data.isalnum():
        return False
    return True
