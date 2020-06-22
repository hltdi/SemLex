"""
Creating cross-linguistic lexical resources for ES languages.
"""

import hm
import re
import random

FS = hm.morpho.fst.FeatStruct

##AMgen = None
##TIgen = None
##KSgen = None
##TEgen = None

AMTIfeats = "[sb=[-fem]]"
KSfeats = "[sg=m]"

TI = None
TIinit = "[cj1=None,cj2=None,-d,-neg,pos=v,pp=None,-rel,sb=[-p1,-p2,-plr],-sub,tm=prf,-yn]"
TE = None
TEinit = "[-neg,op=None,pos=v,-rel,sn=1,sp=3,-sub,tm=prf]"
KS = None

R2C = {}

ENG = re.compile("([a-zA-Z,\-() ]+)")
PAREN = re.compile(r"(\(.+?\))")

# root translation files
ENTRY_SEP = '%'
LANG_SEP = '\n'
ITEM_SEP = ';'
LANG_ITEM_SEP = ':'
ELEMENT_SEP = '.'
RC_SEP = "|"
COMMENT = '#'
# separates voice and aspect features
VA_SEP = ","

GEN = {"am": None, "ti": None, "te": None, "ks": None}
LONG = {"am": "am", "ti": "ti", "te": "tig", "ks": "gru"}

LEXFEAT = {"[as=None,vc=[-ps,-cs]]": "",
           "[as=smp,vc=smp]": "",
           "[as=smp,vc=ps]": "ps",
           "[as=None,vc=[+ps,-cs]]": "ps",
           "[as=smp,vc=tr]": "tr",
           "[as=None,vc=[-ps,+cs]": "tr",
           "[as=smp,vc=cs]": "cs",
           "[as=None,vc=[+ps,+cs]": "cs",
           "[as=it,vc=smp]": ",it",
           "[as=it,vc=[-ps,-cs]]": ",it",
           "[as=it,vc=ps]": "ps,it",
           "[as=it,vc=[+ps,-cs]": "ps,it",
           "[as=it,vc=cs]": "cs,it",
           "[as=it,vc=[-ps,+cs]]": "cs,it",
           # this doesn't occur yet; it's basically cls=C
           "[as=rc,vc=smp]": ",rc",
           "[as=rc,vc=[-ps,-cs]]": ",rc",
           "[as=rc,vc=ps]": "ps,rc",
           "[as=rc,vc=[+ps,-cs]]": "ps,rc",
           "[as=rc,vc=tr]": "tr,rc",
           # this doesn't occur yet; mismatch between AmTi and KsTeCh
           "[as=rc,vc=[-ps,+cs]]": "tr,rc",
           # this doesn't occur yet, though it should
           "[as=rc,vc=cs]": "cs,rc",
           "[as=rc,vc=[+ps,+cs]]": "cs,rc"
           }

At_LEXFEAT = {'': "[as=smp,vc=smp]",
              'imp': "[as=smp,vc=smp,ob=[+xpl]]",
              'trimp': "[as=smp,vc=tr,ob=[+xpl]]",
              'psimp': "[as=smp,vc=ps,ob=[+xpl]]",
              'csimp': "[as=smp,vc=cs,ob=[+xpl]]",
              'ps': "[as=smp,vc=ps]",
              'tr': "[as=smp,vc=tr]",
              'cs': "[as=smp,vc=cs]",
              'it': "[as=it,vc=smp]",
              'trit': "[as=it,vc=tr]",
              'psit': "[as=it,vc=ps]",
              'psitimp': "[as=it,vc=ps,ob=[+xpl]]",
              'csit': "[as=it,vc=cs]",
              'psrc': "[as=rc,vc=ps]",
              'rctr': "[as=rc,vc=tr]",
              '3sf': "[sb=[+fem]]",
              '3sm': "[sb=[-fem]]",
              '1s': "[sb=[+p1,-p2,-plr]]",
              '1p': "[sb=[+p1,-p2,+plr]]",
              '3pf': "[sb=[+fem,-p1,-p2,+plr]]",
              '3pm': "[sb=[-fem,-p1,-p2,+plr]]",
              '2sf': "[sb=[+fem,-p1,+p2,-plr]]",
              '2sm': "[sb=[-fem,-p1,+p2,-plr]]",
              '2pf': "[sb=[+fem,-p1,+p2,+plr]]",
              '2pm': "[sb=[-fem,-p1,+p2,+plr]]",
              '0o': "[ob=[-expl,-xpl]]",
              '3smo': "[ob=[+expl,-fem,-p1,-p2,-plr,-prp,+xpl]]",
              '1so': "[ob=[+expl,+p1,-p2,-plr,-prp,+xpl]]",
              '1po': "[ob=[+expl,+p1,-p2,+plr,-prp,+xpl]]",
              '3sfo': "[ob=[+expl,+fem,-p1,-p2,-plr,-prp,+xpl]]",
              '2pmo': "[ob=[+expl,-fem,-p1,+p2,+plr,-prp,+xpl]]",
              '3pfo': "[ob=[+expl,+fem,-p1,-p2,+plr,-prp,+xpl]]",
              'prf': "[tm=prf]",
              'imf': "[tm=imf]",
              'j_i': "[tm=j_i]",
              'neg': "[+neg]",
              'aff': "[-neg]",
              'rel': "[+rel,+sub]",
              'main': "[-rel,-sub]"
            }

KT_LEXFEAT = {"": "[as=None,vc=[-ps,-cs]]",
              "imp": "[as=None,vc=[-ps,-cs],sp=3,sg=m]",
              "ps": "[as=None,vc=[+ps,-cs]]",
              "psimp": "[as=None,vc=[+ps,-cs],sp=3,sg=m]",
              "psit": "[as=it,vc=[+ps,-cs]]",
              "psitimp": "[as=it,vc=[+ps,-cs],sp=3,sg=m]",
              "psrc": "[as=rc,vc=[+ps,-cs]]",
              "it": "[as=it,vc=[-ps,-cs]]",
              "cs": "[as=None,vc=[-ps,+cs]]",
              "csimp": "[as=None,vc=[-ps,+cs],sp=3,sg=m]",
              "csit": "[as=it,vc=[-ps,+cs]]",
              "csps": "[as=None,vc=[+ps,+cs]]",
              "cspsimp": "[as=None,vc=[+ps,+cs],sp=3,sg=m]",
              "cspsit": "[as=it,vc=[+ps,+cs]]",
              "cspsrc": "[as=rc,vc=[+ps,+cs]]",
              '3sf': "[sg=f]",
              '3sm': "[sg=m]",
              '1s': "[sp=1,sn=1]",
              '1p': "[sp=1,sn=2]",
              '3pf': "[sp=3,sn=2,sg=f]",
              '3pm': "[sp=3,sn=2,sg=m]",
              '2sf': "[sp=2,sn=1,sg=f]",
              '2sm': "[sp=2,sn=1,sg=m]",
              '2pf': "[sp=2,sn=2,sg=f]",
              '2pm': "[sp=2,sn=2,sg=m]",
              '0o': "[op=None]",
              '3smo': "[op=3,on=1,og=m]",
              '1so': "[op=1,on=1]",
              '1po': "[op=1,on=2]",
              '3sfo': "[op=3,on=1,og=f]",
              '2pmo': "[op=2,on=2,og=m]",
              '3pfo': "[op=3,on=2,og=f]",
              'prf': "[tm=prf]",
              'imf': "[tm=imf]",
              'j_i': "[tm=j_i]",
              'neg': "[+neg]",
              'aff': "[-neg]",
              'rel': "[+rel,+sub]",
              'main': "[-rel,-sub]"             
              }

# feature choices for different valence categories
# car of each value is categories features, cadr is a list of sets of features to avoid
VAL_FEATS = \
 {
#         TAM            SBJ                    OBJ                  POL             SUB
  0: [ [ ('prf', 'imf'), ('3sm',),             ('0o',),             ('neg', 'aff'), ('main', 'rel')],
       [] ],
  1: [ [ ('prf', 'imf'), ('3sm',),             ('3smo', '3sfo'),    ('neg', 'aff'), ('main', 'rel')],
       [] ],
  2: [ [ ('prf', 'imf'), ('1s', '3sm', '2pf'), ('0o',),             ('neg', 'aff'), ('main', 'rel')],
       [] ],
  3: [ [ ('prf', 'imf'), ('1s', '3sm', '2pf'), ('3smo', '3sfo'),    ('neg', 'aff'), ('main', 'rel')],
       [] ]
  }

CATS = {'tam': ['imf', 'j_i', 'prf'],
        'sbj': ['1s', '3sf', '2pm'],
        'obj': ['3smo'],
        'pol': ['neg', 'aff'],
        'sub': ['main', 'rel']}

geezify = None

### Generating data

### Valency categories for root|class-asp/vc combinations.
### 0: intransitive: only 3s subjects (only 3sm for 0 valency verbs like <znb|A>)
### 1: impersonal: 3sm subject, any object
### 2: intransitive: all subjects
### 3: transitive: all subjects, only 3s objects for some

def test(vc=3, minimum=4):
    gen_vc_simroots(vc, minimum=minimum)

def convert_data_file(vc, minimum=4, lang_list='ti,am,te,ks', file_comment=''):
    filename = "data_{}_{}".format(vc, minimum)
    languages = lang_list.split(',')
    with open(filename + "_mt.txt", 'w', encoding='utf8') as outf:
        if file_comment:
            print(file_comment, file=outf)
        print("languages = {}".format(lang_list), file=outf)
        print(file=outf)
        with open(filename + ".txt", encoding='utf8') as inf:
            contents = inf.read().split('## ')
            for group in contents:
                if not group:
                    continue
                lines = group.split('\n')
                comment = lines[0]
                print("## {}".format(comment), file=outf)
                forms = [l.split(':') for l in lines[1:] if l]
                forms = dict(forms)
                for lang in languages:
                    lform = forms.get(lang)
                    if not lform:
                        print(file=outf)
                    else:
                        print(lform, file=outf)

def gen_featcombs(vc):
    """Return all feature combinations for valency category vc."""
    featlist, avoidfeats = VAL_FEATS[vc]
    combs = allcombs(featlist, avoidfeats)
    return combs

def gen_vc_simroots(vc, simroots=None, minimum=4, write=True):
    """Given a dict of simroots and a valency category,
    find (or generate) a set of features to use, and
    generate the output word."""
    if not simroots:
        simroots = read_vcats_simroots(minimum)
    featcombs = gen_featcombs(vc)
    forms = []
    for valcat, roots in simroots:
        if valcat == vc:
            for fc in featcombs:
                fc_string = ','.join(fc)
                transset = [fc_string]
                for lang, root in roots.items():
                    rc, f = root.split('.')
                    if len(transset) == 1:
                        # Languages share common asp/vc feature
                        transset.append(f)
                    r, c = rc.split('|')
                    form = gen(lang, r, f, c, morefeats=fc)
                    if form:
                        l, rom, gz = form
                        transset.append((l, rom, gz, rc))
                forms.append(transset)
    if write:
        with open("data_{}_{}.txt".format(vc, minimum), 'w', encoding='utf8') as file:
            for transset in forms:
                f0 = transset[0]
                f1 = transset[1]
                print("## {};{}".format(f1, f0), file=file)
                for l, rom, g, root in transset[2:]:
                    print("{}:{} # {};{}".format(l, rom, root, g), file=file)
    else:
        return forms

def read_vcats_simroots(minimum=1):
    result = []
    with open("eatTk_match2.txt", encoding='utf8') as file:
        for line in file:
            val, roots = line.split(' ; ')
            val = int(val)
            roots = eval(roots)
            if len(roots) < minimum:
                continue
            result.append((val, roots))
    return result

def random_val_feats(valency):
    """Assign random values to different feature categories."""
    feats = VAL_FEATS[valency]
    return [random.choice(f) for f in feats]

def read_vcats():
    """Read in Am root-asp/vc categories."""
    cats = {}
    with open("am_vcats.txt", encoding='utf8') as file:
        for line in file:
            line = line.strip()
            rc, fc = line.split(';')
            cats[rc] = eval(fc)
    return cats

def cat_similar(write="eatTk_match2.txt"):
    cats = read_vcats()
    sim = read_similar()
    new = []
    uncat = []
    for s in sim:
        # s a lang:root|cls.feats dict
        if len(s) < 3:
            continue
        if not 'am' in s:
#            print("No cat for {}".format(s))
            uncat.append(s)
            continue
        am_rcf = s['am']
        am_rc, am_f = am_rcf.split('.')
        a_cats = cats.get(am_rc)
        if not a_cats:
            print("No Am cat for {}_{}".format(am_rc, am_f))
            uncat.append(s)
            continue
        am_f = am_f.replace(',', '')
        c = -1
        if am_f == 'imp':
            c = 1
        else:
            c = a_cats.get(am_f, None)
            if c is None:
                print("No cat for {} in {} ({})".format(am_f, a_cats, am_rcf))
        new.append((c, s))
    return new, uncat

def convert_feat(lang, feat):
    """Convert something like '1s' to something like '[sp=1,sn=1]'."""
    if lang in ('am', 'ti'):
        return At_LEXFEAT.get(feat)
    else:
        return KT_LEXFEAT.get(feat)

def read_similar(file="eatTk_match.txt"):
    res = []
    with open(file, encoding='utf8') as f:
        for line in f:
            langs = line.strip().split(';')
            langs = [l.split(':') for l in langs]
            langs = dict(langs)
            res.append(langs)
    return res

def gen_all(ldict, feats):
    """Given a language: root/features dict (from read_similar())
    and a set of additional features, generate all of the output
    forms.
    """
    res = []
    for lang, rootfeats in ldict.items():
        root, f = rootfeats.split('.')
        r, c = root.split('|')
        out = gen(lang, r, f, c, morefeats=feats)
        res.append(out)
    return res

def gen(lang, root, feats, cls='A', morefeats=None):
    load(lang)
    if morefeats:
        if not isinstance(morefeats, list):
            morefeats = [morefeats]
        morefeats = [convert_feat(lang, mf) for mf in morefeats]
    feats = get_feats(lang, feats, cls, morefeats)
    output = generate(lang, root, feats)
    if output:
        rom = output[0][0]
        geez = geezify(rom, LONG[lang], gemination=True)
        return lang, rom, geez
    else:
        print("Couldn't generate {} {}, {}".format(lang, root, feats.__repr__()))

def get_generator(lang):
    return GEN[lang]

def set_generator(lang, gen):
    GEN[lang] = gen

def generate(lang, root, feats):
    generator = get_generator(lang)
    if not generator:
        return
    return generator.gen(root, update_feats=feats, phon=True if lang in ('am', 'ti') else False)

def load(language):
    global geezify
    if language == 'am' and not get_generator('am'):
        set_generator('am', hm.morpho.get_language('amh', segment=False, phon=True).morphology['v'])
    elif language == 'ti' and not get_generator('ti'):
        set_generator('ti', hm.morpho.get_language('ti', phon=True).morphology['v'])
    elif language == 'ks' and not get_generator('ks'):
        set_generator('ks', hm.morpho.get_language('gru').morphology['v'])
    elif language == 'te' and not get_generator('te'):
        set_generator('te', hm.morpho.get_language('tig').morphology['v'])
    if not geezify:
        geezify = hm.morpho.geez.geezify

# Utilities

def allcombs(seqs, avoid=None):
    """Returns a list of all sequences consisting of one element from each of seqs.
    This could also be done with itertools.product.
    avoid is None or a list of sets of elements to be avoided.
    [From HornMorpho/.../utils.py.]
    """
    def avoid_comb(comb):
        # See if any of the sets in avoid are in comb.
        if not avoid:
            return False
        for a in avoid:
            if a.issubset(comb):
                return True
        return False
    if not seqs:
        return []
    res = [[x] for x in seqs[0]]
    for item in seqs[1:]:
        for i in range(len(res)-1, -1, -1):
            rr = res[i]
            res[i:i+1] = [(rr + [itemitem]) for itemitem in item if not avoid_comb((rr + [itemitem]))]
    return res

def get_feats(lang, feats, cls='A', morefeats=None):
    basic = FS("[cls={}]".format(cls))
    feats = convert_feat(lang, feats)
#    feats = At_LEXFEAT.get(feats) if lang in ('am', 'ti') else KT_LEXFEAT.get(feats)
    feats = FS(feats)
    basic.update(feats)
    if morefeats:
        for mf in morefeats:
            mf = FS(mf)
            basic.update(mf)
    return basic

def get_new_roots():
    roots = {'ks': get_ks_roots(),
             'am': get_am_roots(),
             'ti': get_ti_roots(),
             'te': get_te_roots()}
    new_roots = {'am': [], 'ti': [], 'ks': [], 'te': []}
    with open("eatTk_match.txt", encoding='utf8') as file:
        for line in file:
            for l in line.split(';'):
                language, rcf = l.split(':')
                rc, f = rcf.split('.')
                r, c = rc.split('|')
                if (r, c) not in roots[language]:
                    if (r, c) not in new_roots[language]:
                        new_roots[language].append((r, c))
    return new_roots

def get_ti_roots():
    roots = []
    with open("../HornMorpho/hm/languages/ti/lex/v_root.lex") as file:
        for line in file:
            line = line.split('#')[0]
            if not line:
                continue
            line = line.split()
            root = line[0]
            if len(line) == 1 or 'cls' not in line[-1]:
                cls = 'A'
            else:
                cls = line[-1].split('cls=')[1]
                if ',' in cls:
                    cls = cls.split(',')[0]
                else:
                    cls = cls.split(']')[0]
            roots.append((root, cls))
    return roots

def get_am_roots():
    roots = []
    with open("../HornMorpho/hm/languages/amh/lex/v_root.lex") as file:
        for line in file:
            line = line.split('#')[0]
            if not line:
                continue
            line = line.split()
            root = line[0]
            if len(line) == 1 or 'cls' not in line[-1]:
                cls = 'A'
            else:
                cls = line[-1].split('cls=')[1]
                if ',' in cls:
                    cls = cls.split(',')[0]
                else:
                    cls = cls.split(']')[0]
            roots.append((root, cls))
    return roots

def get_te_roots():
    roots = []
    with open("../HornMorpho/hm/languages/tig/lex/v_root.lex") as file:
        for line in file:
            line = line.split('#')[0]
            if not line:
                continue
            line = line.split()
            root = line[0]
            if len(line) == 1 or 'cls' not in line[-1]:
                cls = 'A'
            else:
                cls = line[-1].split('cls=')[1]
                if ',' in cls:
                    cls = cls.split(',')[0]
                else:
                    cls = cls.split(']')[0]
            roots.append((root, cls))
    return roots    

def get_ks_roots():
    roots = []
    with open("../HornMorpho/hm/languages/gru/lex/v_root.lex") as file:
        for line in file:
            line = line.split('#')[0]
            if not line:
                continue
            line = line.split()
            root = line[0]
            if len(line) == 1 or 'cls' not in line[-1]:
                cls = 'A'
            else:
                cls = line[-1].split('cls=')[1]
                if ',' in cls:
                    cls = cls.split(',')[0]
                else:
                    cls = cls.split(']')[0]
            roots.append((root, cls))
    return roots    

def find_similar(entries, write="eatTk_match.txt"):
    result = []
    for entry in entries:
        rcs = []
        matches = {}
        for language, forms in entry.items():
            rcs1 = [language, []]
            if language == 'en':
                continue
            # get rid of comments
            forms = forms.split(COMMENT)[0]
            for form in forms.split(';'):
                rc, feats, geez = form.split('.')
                root, cls = rc.split('|')
                rcs1[1].append((root, cls, feats))
            rcs.append(rcs1)
        matched = []
        for index, (language, rc) in enumerate(rcs[:-1]):
            for rrcc in rc:
                matched1 = False
                if (language, rrcc) in matched:
                    continue
                for language1, rc1 in rcs[index+1:]:
                    for rrcc1 in rc1:
                        if (language1, rrcc1) in matched:
                            continue
                        if root_match(language, rrcc, language1, rrcc1):
                            matched.append((language1, rrcc1))
                            matched1 = True
                            if (language, rrcc) in matches:
                                matches[(language, rrcc)].append((language1, rrcc1))
                            else:
                                matches[(language, rrcc)] = [(language1, rrcc1)]
                if matched1:
                    matched.append((language, rrcc))
        if matches:
            # Convert the dict into a list of lists
            matches = list(matches.items())
            matches = [m[1] + [m[0]] for m in matches]
##            # There might already be a longer set including these roots
##            for matches1 in matches:
##                l1 = matches1[0][0]
##                r1 = matches1[0][1][0]
##                ignore = False
##                for res in result:
##                    for res1 in res:
##                        if res1[0] == l1 and res1[1][0] == r1:
##                            if len(res) >= len(matches1):
##                                ignore = True
##                                break
##                    if ignore:
##                        break
##                if not ignore:
##                    result.append(matches1)
            result.extend(matches)
    # Sort the results by length
    result.sort(key=lambda r: len(r), reverse=True)
    if write:
        with open(write, 'w', encoding='utf8') as file:
            for result1 in result:
                result1 = ["{}:{}|{}.{}".format(r[0], r[1][0], r[1][1], r[1][2]) for r in result1]
                print(";".join(result1), file=file)
    else:
        return result

def feat_match(f1, f2):
    if f1 == f2:
        return True
    if 'cs' in f1:
        f1 = f1.replace('cs', 'tr')
        if f1 == f2:
            return True
    elif 'cs' in f2:
        f2 = f2.replace('cs', 'tr')
        if f1 == f2:
            return True
    return False

def root_match(l1, rc1, l2, rc2):
#    print("Checking {}:{} and {}:{}".format(l1, rc1, l2, rc2))
    def root_match1(root1, cls1, feat1, root2, cls2, feat2):
        return root1 == root2 and cls1 == cls2 and feat_match(feat1, feat2)
    r1, c1, f1 = rc1
    r2, c2, f2 = rc2
    if l1 == 'am':
        r2 = r2.replace("`", "'").replace("h", "'").replace("H", "'")
    if l2 == 'am':
        r1 = r1.replace("`", "'").replace("h", "'").replace("H", "'")
    if root_match1(r1, c1, f1, r2, c2, f2):
        return True
    if 'S' in r1:
        r1 = r1.replace('S', 'T')
        if root_match1(r1, c1, f1, r2, c2, f2):
            return True
    if 'S' in r2:
        r2 = r2.replace('S', 'T')
        if root_match1(r1, c1, f1, r2, c2, f2):
            return True
    if 'x' in r1:
        r1 = r1.replace('x', 's')
        if root_match1(r1, c1, f1, r2, c2, f2):
            return True
    if 'x' in r2:
        r2 = r2.replace('x', 's')
        if root_match1(r1, c1, f1, r2, c2, f2):
            return True
    if 'W' in r1 and 'W' in c2:
        r1 = r1.replace('W', '')
        if root_match1(r1, c1, f1, r2, c2, f2):
            return True
    if 'W' in r2 and 'W' in c1:
        r2 = r2.replace('W', '')
        if root_match1(r1, c2, f1, r2, c2, f2):
            return True

def write(entries, path="eatTk4.txt", sort=True):
    """
    Write a lexicon, a list of dict entries, to file.
    """
    with open(path, 'w', encoding='utf8') as file:
        if sort:
            # Alphabetize by English
            entries.sort(key=lambda e: e['en'])
        for entry in entries:
            print(ENTRY_SEP, file=file)
            lgs = []
            for lg, items in entry.items():
                lgs.append("{}{}{}".format(lg, LANG_ITEM_SEP, items))
            lgs = LANG_SEP.join(lgs)
            print(lgs, file=file)

def read(filename='eatTk4.txt', minimum=0, maximum=0, include='',
         update=None):
    """
    Read a lexicon in from a file as a list of dicts.
    minimum and maximum constrain the number of languages (other than English) that
    must be in an entry.
    include is a language abbreviation string ('am', 'ti', 'te', 'ks') specifying
    a language that must be in an entry.
    """
    with open(filename, encoding='utf8') as file:
        contents = file.read()
        return [e for e in read_entries(contents, minimum=minimum, maximum=maximum, include=include, update=update) if e]

def read_entries(entries, minimum=0, maximum=0, include='', update=None):
    e = []
    entries = entries.split(ENTRY_SEP)
    return [read_entry(entry.strip(), minimum=minimum, maximum=maximum, include=include, update=update) for entry in entries if entry]

def read_entry(entry, minimum=0, maximum=0, include='', update=None):
    edict = {}
    langs = entry.split(LANG_SEP)
    langs = [lang.strip() for lang in langs]
#    print("Read entry: {}, include: {}".format(entry, include))
    # number of languages, excluding English
    # current maximum: 4 (am, ti, te, ks)
    nlangs = len(langs) - 1
    if maximum and nlangs > maximum:
        return
    if minimum and nlangs < minimum:
        return
    if include and not any([lang.startswith(include) for lang in langs]):
        return
    for lang in langs:
#        print(lang)
        # separate off possible commments
        item, comment = sep_comment(lang)
        item_split = item.split(LANG_ITEM_SEP)
        if len(item_split) != 2:
            print("Something wrong with {} # {}".format(item, comment))
        l, items = item_split
        if update and l in update:
            # Update the items for this language
            new_items = []
            for item in items.split(ITEM_SEP):
                elements = item.split(ELEMENT_SEP)
                root = elements[0]
                if l in ['ks', 'te']:
                    features = elements[-2]
                    features = features.split('cls=')
                    if len(features) != 2:
                        print(root, features)
                    features = features[1].split(',')
                    cls = features[0]
                    features = '[' + ','.join(features[1:])
                    item = "{}{}{}{}{}{}{}".format(root, RC_SEP, cls, ELEMENT_SEP, features, ELEMENT_SEP, elements[-1])
                else:
                    new_root, cls = convert_root(root, lg=l)
                    rc = "{}{}{}".format(new_root, RC_SEP, cls)
                    item = "{}{}{}".format(rc, ELEMENT_SEP, ELEMENT_SEP.join(elements[1:]))
                if item not in new_items:
                    new_items.append(item)
            items = ITEM_SEP.join(new_items)
        if comment:
            items = "{}{}{}".format(items, COMMENT, comment)
        edict[l] = items
    return edict

def sep_comment(item):
    """Separate an item into the content and a possible comment."""
    item_comment = item.split(COMMENT)
    comment = ''
    if len(item_comment) > 1:
        # Length should be 2 but could be greater if there are multiple comment characters
        item = item_comment[0].strip()
        comment = item_comment[1].strip()
    return item, comment

def roots2classes(lg='ti'):
    if lg in R2C:
        return R2C[lg]
    rc = {}
    path = ''
    if lg == 'ti':
        path = "../LingData/Ti/roots2class.txt"
    elif lg == 'am':
        path = "../LingData/Am/roots2class.txt"
    if path:
        with open(path, encoding='utf8') as file:
            for line in file:
                root, classes = line.strip().split()
                rc[root] = classes.split(',')
        R2C[lg] = rc
        return rc
    else:
        print("Don't know how to update {}".format(lg))

def convert_root(root, rc=None, lg='ti'):
    """Convert an old-style Am or Ti root string to a new root, class pair."""
    rc = rc or roots2classes(lg=lg)
    reduced = root.replace('_', '').replace('|', '').replace('a', '')
    entry = rc.get(reduced)
    if not entry:
        print("Something wrong: {} not in lexicon".format(root))
        return None, ''
    elif len(entry) == 1:
        return reduced, entry[0]
    else:
        if '_' in root and 'B' in entry:
            return reduced, 'B'
        elif 'a' in root:
            if 'C' in entry:
                return reduced, 'C'
            elif 'F' in entry:
                return reduced, 'F'
            elif 'J' in entry:
                return reduced, 'J'
        elif 'A' in entry:
            return reduced, 'A'
        elif 'E' in entry:
            return reduced, 'E'
        else:
            print("Something wrong with {}, {}".format(root, entry))
            return None, ''

def elim_paren(string):
    return PAREN.sub("", string)

def proc_en_te_ti():
    results = []
    with open("Tet1.txt", encoding='utf8') as file:
        for line in file:
            te, en, ti = line.split(';;')
            ti = ti.strip()
            tigeez, tirf = ti.split(';')
            tir, tif = tirf.split(":")
            tigeezgem = ti_gen(tir, tif, phon=True)
            if not tigeezgem:
                print("Couldn't generate {}".format(tirf))
            tigeezgem = tigeezgem[1]
            ti = "{}:{}".format(tigeezgem, tirf)
            results.append((en, te, ti))
    results.sort()
    with open("eTt1.txt", 'w', encoding='utf8') as file:
        for en, te, ti in results:
            print(ENTRY_SEP, file=file)
            print("en;{};;\nte;{};;\nti;{}".format(en, te, ti), file=file)
    return results

def proc_te_ti():
    load_anal('ti')
    load_anal('te')
#    estems = eng_verbs()
    epp = eng_pp()
    result = []
#    new_te = []
    with open("te_tr.txt", encoding='utf8') as file:
        for line in file:
            line = line.split("\t")
            if len(line) != 3:
                print("{} -- {}".format(len(line), line.strip()))
            te, ti, en = line
            en = en.strip()
            # analyze Ti word
            ti = ti.split("፡")
            en = en.split(',')
            enproc = []
            # Tigre
            B = False
            te = te.strip()
            if "'" in te:
                B = True
            te = te.replace("'", "")
            teanal = te_anal(te, guess=False)
            if not teanal:
                print("Couldn't analyze {}".format(te))
            if len(teanal) > 1:
                print("Multiple analyses for {}: {}".format(te, teanal))
            teanal = teanal[0]
#                teanal = te_anal(te, guess=True)
#                new_te.append((te, teanal))
#                if not teanal:
#                    print("Couldn't analyze {}".format(te))
            # English
            for ee in en:
                ee = ee.strip()
                eesplit = ee.split()
#                if eesplit[0] in estems:
                root = eesplit[0] + '_v'
                rest = eesplit[1:]
                enproc.append(' '.join([root] + rest))
                if rest and root == 'be_v':
                    next_word = rest[0].strip()
                    if next_word in epp:
                        enprocps = epp[next_word] + "_v[ps]"
                        enprocps = ' '.join([enprocps] + rest[1:])
#                        print("Creating passive entry {}".format(enprocps))
                        enproc.append(enprocps)
#                else:
#                    print("{} not an English verb".format(eesplit[0]))
            # Tigrinya
            for tti in ti:
                B = False
                tti = tti.strip()
                if tti and not ' ' in tti:
                    if "'" in tti:
                        B = True
                    tti = tti.replace("'", "")
                    tianal = ti_anal(tti, guess=False)
                    if len(tianal) > 1:
#                        print("Multiple {}".format(tianal))
                        if B:
                            tianal = [a for a in tianal if '_' in a]
                        else:
                            tianal = [a for a in tianal if '_' not in a]
                    tianal = tianal[0]
#                    tigem = ti_gen(tianal[0], tianal[1], True)
#                    if not tigem:
#                        print("Couldn't generate {}".format(tianal)
                    for enp in enproc:
                        res = "{};{};;{};;{};{}".format(te, teanal, enp, tti, tianal)
                        result.append(res)
    result.sort()
    with open("Tet1.txt", 'w', encoding='utf8') as file:
        for r in result:
            print(r, file=file)
    return result

def update_ti_roots():
    new_roots = []
    old_roots = get_ti_roots()
    with open("tmp.txt", encoding='utf8') as file:
        for line in file:
            g, e, rf = line.split(';')
            r, f = rf.split(':')
            r = r.strip()
            if r in old_roots:
                print("{} already in lexicon".format(r))
            else:
                new_roots.append(r)
    return new_roots

def get_new_ti_roots():
    old = get_ti_roots()
    roots = []
    with open("ti_new.txt", encoding='utf8') as file:
        for line in file:
            ti, en = line.split(';;')
            tigeez, tirf = ti.split(';')
            tiroot = tirf.split(':')[0]
            if tiroot in old:
                print("{} is already in lexicon".format(tiroot))
            else:
                roots.append(tiroot)
    return roots

def load_anal(language):
    global TI
    global TE
    if language == 'ti' and not TI:
        TI = hm.morpho.get_language('ti', segment=False)
    elif language == 'te' and not TE:
        TE = hm.morpho.get_language('tig')

def ti_anal(word, guess=False):
    load_anal('ti')
    anal = TI.anal_word(word, init_weight=TIinit, preproc=True, guess=guess)
    if anal:
        return [tiam_anal2string(a) for a in anal]
#        return [(x[0], x[1].get('as'), x[1].get('vc')) for x in anal]

def te_anal(word, guess=False):
    load_anal('te')
    anal = TE.anal_word(word, init_weight=TEinit, preproc=True, guess=guess)
    if anal:
        return [teks_anal2string(a) for a in anal]

def teks_anal2string(anal):
    root = anal[0]
    feats = anal[1]
    ob = feats.get('op')
    featstring = "[cls={},as={},vc={}".format(feats.get('cls'), feats.get('as', 'smp').__repr__(), feats.get('vc', 'smp').__repr__())
    if ob:
        featstring += ",op=3,on=1,og=m"
    featstring += "]"
    return "{}:{}".format(root, featstring)

def tiam_anal2string(anal):
    root = anal[0]
    feats = anal[1]
    ob = feats['ob']['xpl']
    featstring = "[as={},vc={}".format(feats.get('as', 'smp'), feats.get('vc', 'smp'))
    if ob:
        featstring += ",ob=[+xpl]"
    featstring += "]"
    return "{}:{}".format(root, featstring)

def eng_verbs():
    verbs = []
    with open("../LingData/En/v_analyzed.lex") as file:
        for line in file:
            if "tam=stm" in line:
                stem = line.split()[0]
                verbs.append(stem)
    verbs.sort()
    return verbs

def eng_pp():
    pp = {}
    with open("../LingData/En/v_analyzed.lex") as file:
        for line in file:
            if "tam=pp" in line:
                forms = line.split()
                pp[forms[0]] = forms[1]
    return pp


def eng_passive(defn, pps=None):
    pps = pps or eng_pp()
    dsplit = defn.split()
    if dsplit[0] == "be_v":
        if len(dsplit) > 1 and dsplit[1] in pps:
            ppstem = pps[dsplit[1]]
            return ppstem + "_v[ps]"
    return ''
