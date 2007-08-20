## TLS Motion Determination (TLSMD)
## Copyright 2002-2005 by TLSMD Development Group (see AUTHORS file)
## This code is part of the TLSMD distribution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.

import os
import sys
import subprocess
import traceback

from django.conf import settings
from django.core.mail import send_mail, mail_admins


def SendEmail(address, subject, body):
    send_mail(subject, body, settings.SERVER_EMAIL, [address], fail_silently=False)


def SendTracebackEmail(context):
    mail_admins(context, traceback.format_exc())
