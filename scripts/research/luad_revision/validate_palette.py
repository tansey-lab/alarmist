"""Python port of the dataviz skill's validate_palette.js six checks.

Needed because node on this machine is v10 (2018) and cannot load the shipped ES module.
Thresholds, the Machado-Oliveira-Fernandes (2009) severity-1.0 CVD matrices and the OKLab
conversion are copied verbatim from scripts/validate_palette.js.
"""
import math, sys

BAND = {'light': (0.43, 0.77), 'dark': (0.48, 0.67)}
CHROMA_FLOOR = 0.10
CVD_TARGET, CVD_FLOOR = 8.0, 6.0
NORMAL_FLOOR = 15.0
CONTRAST_MIN = 3.0
DEFAULT_SURFACE = {'light': '#fcfcfb', 'dark': '#1a1a19'}
MACHADO = {
 'protan': [[0.152286,1.052583,-0.204868],[0.114503,0.786281,0.099216],[-0.003882,-0.048116,1.051998]],
 'deutan': [[0.367322,0.860646,-0.227968],[0.280085,0.672501,0.047413],[-0.011820,0.042940,0.968881]],
 'tritan': [[1.255528,-0.076749,-0.178779],[-0.078411,0.930809,0.147602],[0.004733,0.691367,0.303900]],
}
def hex2srgb(h):
    h = h.strip().lstrip('#')
    return [int(h[i:i+2], 16)/255 for i in (0,2,4)]
def s2lin(c): return c/12.92 if c <= 0.04045 else ((c+0.055)/1.055)**2.4
def lin(h): return [s2lin(c) for c in hex2srgb(h)]
def relLum(h):
    r,g,b = lin(h); return 0.2126*r + 0.7152*g + 0.0722*b
def contrast(a,b):
    hi,lo = sorted([relLum(a),relLum(b)], reverse=True); return (hi+0.05)/(lo+0.05)
def oklab_from_lin(rgb):
    r,g,b = rgb
    l = (0.4122214708*r + 0.5363325363*g + 0.0514459929*b)**(1/3)
    m = (0.2119034982*r + 0.6806995451*g + 0.1073969566*b)**(1/3)
    s = (0.0883024619*r + 0.2817188376*g + 0.6299787005*b)**(1/3)
    return (0.2104542553*l + 0.7936177850*m - 0.0040720468*s,
            1.9779984951*l - 2.4285922050*m + 0.4505937099*s,
            0.0259040371*l + 0.7827717662*m - 0.8086757660*s)
def oklch(h):
    L,a,b = oklab_from_lin(lin(h)); return L, math.hypot(a,b)
def simulate(h, kind):
    r,g,b = lin(h); M = MACHADO[kind]
    return [min(1,max(0, M[i][0]*r + M[i][1]*g + M[i][2]*b)) for i in range(3)]
def deltaE(h1,h2,kind=None):
    a = oklab_from_lin(simulate(h1,kind) if kind else lin(h1))
    b = oklab_from_lin(simulate(h2,kind) if kind else lin(h2))
    return 100*math.dist(a,b)

pal = [c.strip() for c in sys.argv[1].split(',') if c.strip()]
mode = sys.argv[2] if len(sys.argv) > 2 else 'light'
allpairs = '--pairs' in sys.argv and 'all' in sys.argv
surface = DEFAULT_SURFACE[mode]
fail = False
print(f'palette {pal}   mode {mode}   surface {surface}   '
      f'pairs {"all" if allpairs else "adjacent"}\n')
print(f'{"slot":<10}{"OKLCH L":>9}{"band":>7}{"C":>8}{"chroma":>8}{"contrast":>10}{"":>3}')
lo_b, hi_b = BAND[mode]
for h in pal:
    L, C = oklch(h); ct = contrast(h, surface)
    b_ok = lo_b <= L <= hi_b; c_ok = C >= CHROMA_FLOOR; ct_ok = ct >= CONTRAST_MIN
    fail |= (not b_ok) or (not c_ok)
    print(f'{h:<10}{L:>9.3f}{"PASS" if b_ok else "FAIL":>7}{C:>8.3f}'
          f'{"PASS" if c_ok else "FAIL":>8}{ct:>10.2f}{"" if ct_ok else " WARN":>3}')
pairs = ([(i,j) for i in range(len(pal)) for j in range(i+1,len(pal))] if allpairs
         else [(i,i+1) for i in range(len(pal)-1)])
print(f'\n{"pair":<24}{"normal":>9}{"protan":>9}{"deutan":>9}{"tritan":>9}{"min(p,d)":>10}   verdict')
worst_normal = 1e9
for i,j in pairs:
    n = deltaE(pal[i],pal[j]); p = deltaE(pal[i],pal[j],'protan')
    d = deltaE(pal[i],pal[j],'deutan'); t = deltaE(pal[i],pal[j],'tritan')
    m = min(p,d); worst_normal = min(worst_normal, n)
    v = 'PASS' if m >= CVD_TARGET else ('WARN (floor band)' if m >= CVD_FLOOR else 'FAIL')
    fail |= m < CVD_FLOOR
    print(f'{pal[i]+" / "+pal[j]:<24}{n:>9.1f}{p:>9.1f}{d:>9.1f}{t:>9.1f}{m:>10.1f}   {v}')
nv = 'PASS' if worst_normal >= NORMAL_FLOOR else 'FAIL'
fail |= worst_normal < NORMAL_FLOOR
print(f'\nnormal-vision floor: worst pair {worst_normal:.1f} (needs >= {NORMAL_FLOOR})  {nv}')
print('\nRESULT:', 'FAIL' if fail else 'PASS')
sys.exit(1 if fail else 0)
