"""
Python wrapper around xselect
"""
import subprocess
import time

class Xselect(object):
    def __init__(self, session_name=None):
        self.session_name = 'xsel%d' % int(time.time() % 604800) if session_name is None else session_name
        self.xsel = None

    def __enter__(self):
        args = ['xselect',
                'prefix=%s' % session_name]

        self.xsel = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=None, universal_newlines=True, text=True, bufsize=1)

    def command(self, cmd):
        self.xsel.stdin.write(cmd + '\n')
        self.xsel.stdin.flush()

    def __exit__(self, exc_type, exc_value, traceback):
        self.command('exit')
        self.command('no')
        self.xsel.stdin.close()
        self.xsel.wait()
