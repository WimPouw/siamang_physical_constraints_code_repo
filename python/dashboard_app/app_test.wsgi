import sys
sys.path.insert(0, "/data/www/wimpouw/web/siamang")
sys.path.insert(0, '/data/www/wimpouw/.local/lib/python3.8/site-packages')

from app import app
application = app.server