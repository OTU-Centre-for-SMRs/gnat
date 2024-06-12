#!/usr/bin/env python3

import mms

fs,ss = mms.evaluate('diff(u,t) - div(D * grad(u)) + vec.dot(grad(u)) + alpha * u', 'x * y * t**3', scalars=['D', 'alpha'], vectors=['vec'])
mms.print_fparser(fs)
mms.print_hit(fs, 'force')
mms.print_hit(ss, 'exact')
