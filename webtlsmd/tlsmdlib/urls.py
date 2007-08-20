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

    (r'^tlsmd/static/(?P<path>.*)$',
     'django.views.static.serve',
     {'document_root': os.path.join(settings.TLSMD_ROOT, 'static'), 'show_indexes':True}),
    
    (r'^$', 'tlsmdlib.views.default_redirect'),
)

