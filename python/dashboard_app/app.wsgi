#!/usr/bin/env python3.8

import sys
import logging
import site

site_packages = '/data/home/wimpouw/.local/lib/python3.8/site-packages'
site.addsitedir(site_packages)

sys.path.insert(0, "/data/www/wimpouw/web/siamang")
logging.basicConfig(stream=sys.stderr)

from app import app
application = app.server
