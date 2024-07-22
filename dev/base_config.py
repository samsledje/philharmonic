#!/usr/bin/env python

PROJ_DIR=".."

DEBUG_MODE=True #False #True

def dbg_print(s, *args, **kwargs):
    if DEBUG_MODE or s[:5] != "Flag ":
        print(s, *args, **kwargs)

