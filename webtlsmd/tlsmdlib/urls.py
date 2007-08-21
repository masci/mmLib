import os
from django.conf.urls.defaults import *
from django.conf import settings


urlpatterns = patterns(
    '',

    (r'^tlsmd/$', 'tlsmdlib.views.index'),
    (r'^tlsmd/submit/$', 'tlsmdlib.views.submit'),
    (r'^tlsmd/jobs/$', 'tlsmdlib.views.jobs'),
    (r'^tlsmd/job/(?P<job_id>[A-Za-z0-9_]+)/$', 'tlsmdlib.views.job'),
    (r'^tlsmd/examples/$', 'tlsmdlib.views.examples'),
    (r'^tlsmd/documentation/$', 'tlsmdlib.views.documentation'),

    (r'^accounts/profile/$', 'tlsmdlib.views.my_account'),
    (r'^tlsmd/login/$', 'django.contrib.auth.views.login', {'template_name': 'webtlsmd/login.html'}),
    (r'^tlsmd/logout/$', 'django.contrib.auth.views.logout', {'template_name': 'webtlsmd/logout.html'}),

    (r'^tlsmd/static/(?P<path>.*)$',
     'django.views.static.serve',
     {'document_root': os.path.join(settings.TLSMD_ROOT, 'static')}),

    (r'^tlsmd/jobs/(?P<path>.*)$',
     'django.views.static.serve',
     {'document_root': settings.TLSMD_WORK_DIR, 'show_indexes':True}),
    
    (r'^$', 'tlsmdlib.views.default_redirect'),
)

