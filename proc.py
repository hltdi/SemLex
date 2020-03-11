"""
Creating cross-linguistic lexical resources for ES languages.
"""

import hm
import re

FS = hm.morpho.fst.FeatStruct

AMgen = None
TIgen = None
KSgen = None

AMTIfeats = "[sb=[-fem]]"
#TE = hm.morpho.get_language('te').morphology['v']
#KS = hm.morpho.get_language('ks').morphology['v']
#CH = hm.morpho.get_language('ch').morphology['v']

#AManal = hm.morpho.get_language('am', segment=False)
TI = None
TIinit = "[cj1=None,cj2=None,-d,-neg,pos=v,pp=None,-rel,sb=[-p1,-p2,-plr],-sub,tm=prf,-yn]"
TE = None
TEinit = "[-neg,op=None,pos=v,-rel,sn=1,sp=3,-sub,tm=prf]"
KS = None
KSfeats = "[sg=m]"

ENG = re.compile("([a-zA-Z,\-() ]+)")
PAREN = re.compile(r"(\(.+?\))")

# root translation files
ENTRY_SEP = ';;;'
LANG_SEP = ';;'
ITEM_SEP = ';'
LANG_ITEM_SEP = '::'

AF_CONV = {'': "[as=smp,vc=smp]",
          '3': "[as=smp,vc=smp,ob=[+xpl]]",
          'tr3': "[as=smp,vc=tr,ob=[+xpl]]",
          'ps3': "[as=smp,vc=ps,ob=[+xpl]]",
          'cs3': "[as=smp,vc=cs,ob=[+xpl]]",
          'ps': "[as=smp,vc=ps]",
          'tr': "[as=smp,vc=tr]",
          'cs': "[as=smp,vc=cs]",
          'it': "[as=it,vc=smp]",
          'trit': "[as=it,vc=tr]",
          'psit': "[as=it,vc=ps]",
          'psit3': "[as=it,vc=ps,ob=[+xpl]]",
          'csit': "[as=it,vc=cs]",
          'psrc': "[as=rc,vc=ps]",
          'rctr': "[as=rc,vc=tr]"}

KF_CONV = {"": "[as=None,vc=[-ps,-cs]]",
           "3": "[as=None,vc=[-ps,-cs],sp=3,sg=m]",
           "ps": "[as=None,vc=[+ps,-cs]]",
           "ps3": "[as=None,vc=[+ps,-cs],sp=3,sg=m]",
           "psit": "[as=it,vc=[+ps,-cs]]",
           "psit3": "[as=it,vc=[+ps,-cs],sp=3,sg=m]",
           "psrc": "[as=rc,vc=[+ps,-cs]]",
           "it": "[as=it,vc=[-ps,-cs]]",
           "cs": "[as=None,vc=[-ps,+cs]]",
           "cs3": "[as=None,vc=[-ps,+cs],sp=3,sg=m]",
           "csit": "[as=it,vc=[-ps,+cs]]",
           "csps": "[as=None,vc=[+ps,+cs]]",
           "csps3": "[as=None,vc=[+ps,+cs],sp=3,sg=m]",
           "cspsit": "[as=it,vc=[+ps,+cs]]",
           "cspsrc": "[as=rc,vc=[+ps,+cs]]"}

geezify = None

def merge_ks():
    adict = get_adict()
    ak, ke = get_ak_ek()
#    nokg = []
    emiss = []
    for a, others in adict.items():
        if a in ak:
            k = ak[a]
            if len(others) > 1:
                found = False
                # Multiple Eng translations for Am verb
#                print("{} in AK: {}".format(a, k))
                aengs = list(others.keys())
#                print("  Eng for Am {}".format(list(others.keys())))
                kg = k.split(':')[-1]
                if kg not in ke:
#                    print("  ?? {} not in KE dict".format(kg))
#                    nokg.append((a, k, kg, aengs))
                    pass
                else:
                    kengs = ke[kg]
#                    print("  Eng for Ks {}".format(kengs))
                    for aeng in aengs:
                        for keng in kengs:
                            if keng.startswith(aeng):
#                                print("  a {}, k {}: {}={}".format(a, k, aeng, keng))
                                adict[a][aeng].append(('ks', k))
                                found = True
#                    if not found:
#                        emiss.append((kg, aengs, kengs))
            else:
                eng0 = list(others.keys())[0]
                adict[a][eng0].append(('ks', k))
    return adict

def get_ak_ek():
    ak = {}
    ke = {}
    with open("ak3.txt", encoding='utf8') as file:
        for line in file:
            a, k = line.strip().split('::')
            ak[a] = k
#            kr, ka, kg = k.split(':')
#            if kg in kdict:
#                kdict[kg].append(a)
#            else:
#                kdict[kg] = [a]
    with open("ke1.txt", encoding='utf8') as file:
        for line in file:
            e, k = line.strip().split('::')
            if k in ke:
                ke[k].append(e)
            else:
                ke[k] = [e]
    return ak, ke

##def k2e():
##    kdict = {}
##    kdefs = {}
##    ev = eng_verbs()
##    ep = eng_pp()
##    with open("ak3.txt", encoding='utf8') as file:
##        for line in file:
##            a, k = line.strip().split('::')
##            kr, ka, kg = k.split(':')
##            if kg in kdict:
##                kdict[kg].append(a)
##            else:
##                kdict[kg] = [a]
##    with open("../LingData/Ks/kdict3.txt", encoding='utf8') as file:
##        for line in file:
##            line = line.strip()
##            splitline = line.split()
##            ks = splitline[0]
##            ks = ''.join([x for x in ks if not x.isdigit()])
##            if ks in kdict:
##                eng = ')'.join(line.split(')')[1:])
##                eng = ENG.search(eng)
##                if eng:
##                    eng = eng.group(1).strip()
##                    engs = eng.split(',')
##                    for eng in engs:
##                        eng = elim_paren(eng)
##                        esplit = eng.split()
##                        if esplit:
##                            everb = esplit[0]
##                            if everb in ev:
##                                if ' or ' in eng:
##                                    engdisj = eng.split(' or ')
##                                    # is second word a verb
##                                    engdisj2 = engdisj[1]
##                                    engdisj2split = engdisj2.split()
##                                    everb2 = engdisj2split[0]
##                                    if everb2 in ev:
##                                        # add two verbs
##                                        eng1 = ' '.join([everb + "_v"] + engdisj[0].split()[1:])
##                                        eng2 = ' '.join([everb2 + "_v"] + engdisj2split[1:])
##                                        eng = [eng1, eng2]
##                                    elif everb == 'be':
##                                        # OK, the second disjunct must be an adjective
##                                        eng1 = ' '.join([everb + "_v"] + engdisj[0].split()[1:])
##                                        eng2 = "be_v " + engdisj2
##                                        eng = [eng1, eng2]
##                                        if everb2 in ep:
##                                            # passive
##                                            eng3 = ep[everb2] + "_v[ps]"
##                                            eng.append(eng3)
##                                            edisj1_2 = engdisj[0].split()[1]
##                                            if edisj1_2 in ep:
##                                                eng4 = ep[edisj1_2] + "_v[ps]"
##                                                eng.append(eng4)
##                                    else:
##                                        eng = [eng1]
##                                else:
##                                    eng = [' '.join([everb + "_v"] + esplit[1:])]
##                                    if everb == 'be':
##                                        if esplit[1:][0] in ep:
##                                            eng2 = ep[esplit[1:][0]] + "_v[ps]"
##                                            eng.append(eng2)
##                                if ks in kdefs:
##                                    kdefs[ks].extend(eng)
##                                else:
##                                    kdefs[ks] = eng
##    with open("ke1.txt", 'w', encoding='utf8') as file:
##        for k, eng in kdefs.items():
##            for e in eng:
##                print("{}::{}".format(e, k), file=file)
##    return kdict, kdefs

def get_adict():
    entries = read("eatTk1.txt")
    adict = {}
    for entry in entries:
        aitems = entry.get('am')
        if aitems and len(entry) > 2:
            eitem = entry['en']
            # There has to be a least one entry other than English
            others = [items for items in entry.items() if items[0] != 'am' and items[0] != 'en']
#            others = dict(others)
            aitems = aitems.split(ITEM_SEP)
            for aitem in aitems:
                if aitem in adict:
                    adict[aitem][eitem] = others
                else:
                    adict[aitem] = {eitem: others}
    return adict

##def write_converted_entries(entries, lg, path="eatTk1.txt"):
##    """Write in canonical form (one English item per entry) a dict of entries
##    with keys in another language."""
##    edict = {}
##    for lg_item, eng_items in entries.items():
##        for eng, other_items in eng_items.items():
##            if eng not in edict:
##                edict[eng] = {}
##            eentry = edict[eng]
##            for other_lg, other_item in other_items:
##                eentry[other_lg] = other_item
##            eentry[lg] = lg_item
##    elist = list(edict.items())
##    elist.sort()
##    with open (path, 'w', encoding='utf8') as file:
##        for eng, items in elist:
##            print(ENTRY_SEP, file=file)
##            print("en{}{}{}".format(LANG_ITEM_SEP, eng, LANG_SEP), file=file)
##            others = []
##            for lg, lg_item in items.items():
##                others.append("{}{}{}".format(lg, LANG_ITEM_SEP, lg_item))
##            others = (LANG_SEP + '\n').join(others)
##            print(others, file=file)
##
##def write(entries, path):
##    with open(path, 'w', encoding='utf8') as file:
##        for entry in entries:
##            print(ENTRY_SEP, file=file)
##            lgs = []
##            for lg, items in entry.items():
##                lgs.append("{}{}{}".format(lg, LANG_ITEM_SEP, items))
##            lgs = (LANG_SEP + "\n").join(lgs)
##            print(lgs, file=file)

def read_entry(entry, minimum=0, maximum=0):
    edict = {}
    langs = entry.split(LANG_SEP)
    if maximum and len(langs) > maximum:
        return
    if minimum and len(langs) < minimum:
        return
    for lang in langs:
#        print(lang)
        l, items = lang.strip().split(LANG_ITEM_SEP)
        edict[l] = items
    return edict

def read_entries(entries, minimum=0, maximum=0):
    e = []
    entries = entries.split(ENTRY_SEP)
    return [read_entry(entry.strip(), minimum=minimum, maximum=maximum) for entry in entries if entry]

def read(filename='eatT1.txt', minimum=0, maximum=0):
    with open(filename, encoding='utf8') as file:
        contents = file.read()
        return [e for e in read_entries(contents, minimum=minimum, maximum=maximum) if e]

def elim_paren(string):
    return PAREN.sub("", string)

def proc_ak2():
    entries = {}
    with open("ak2.txt", encoding='utf8') as file:
        for line in file:
            am, ks = line.split('::')
            kroot, kfeats = ks.strip().split('.')
            kcls, kfkey = kfeats.split('_')
            kfeats = KF_CONV[kfkey]
            if not kfeats:
                print("No feats for {}, {}".format(kroot, kfkey))
            kfeats = "[cls={},{}".format(kcls, kfeats[1:])
            kgeez = ks_gen(kroot, kfeats, kcls)
            if kgeez:
                kid = "{}:{}:{}".format(kroot, kfeats, kgeez[-1])
                entries[am] = kid
    entries = list(entries.items())
    entries.sort()
    with open("ak3.txt", 'w', encoding='utf8') as file:
        for key, entry in entries:
            print("{}::{}".format(key, entry), file=file)
#    return entries

##def proc_ak():
##    # Am keyed dict
##    entries = {}
##    with open("ak1.txt", encoding='utf8') as file:
##        for line in file:
##            aroot, trans = line.split(';')
##            trans = eval(trans)
##            for t in trans:
##                afeat, ks = t.split(':')
##                afeat = AF_CONV.get(afeat)
##                if not afeat:
##                    print("No way to convert feats for {}".format(aroot))
##                    return
##                ageez = am_gen(aroot, afeat, True)
##                ageez = ageez[-1]
###                print("aroot {}, afeat {}, ageez {}".format(aroot, afeat, ageez))
##                aid = "{}:{}:{}".format(aroot, afeat, ageez)
##                entries[aid] = ks
##    entries = list(entries.items())
##    entries.sort()
##    with open("ak2.txt", 'w', encoding='utf8') as file:
##        for key, entry in entries:
##            print("{}::{}".format(key, entry), file=file)
###    return entries

##def merge_eatT():
##    # English-keyed dict
##    entries = {}
##    with open("eat1.txt", encoding='utf8') as eat:
##        with open("eTt1.txt", encoding='utf8') as eTt:
##            for entry in eat.read().split(ENTRY_SEP):
##                languages = entry.strip().split(';;')
##                eng = languages[0].split(';')[1]
##                entries[eng] = {}
##                for language in languages[1:]:
##                    forms = language.strip().split(';')
##                    entries[eng][forms[0]] = forms[1:]
##            for entry in eTt.read().split(ENTRY_SEP):
##                languages = entry.strip().split(';;\n')
##                eng = languages[0].split(';')[1]
##                if eng not in entries:
##                    entries[eng] = {}
##                e = entries[eng]
##                for language in languages[1:]:
##                    forms = language.strip().split(';')
##                    l = forms[0]
##                    translation = forms[1:]
##                    if len(translation) == 2:
##                        geez, rf = translation
##                    else:
##                        geez, root, feats = translation[0].split(':')
##                        rf = "{}:{}".format(root, feats)
##                    trans = "{}:{}".format(rf, geez)
##                    if l not in e:
##                        e[l] = [trans]
##                    else:
##                        el = e[l]
##                        if trans not in el:
###                            print("{} already in {}".format(trans, el))
###                        else:
##                            el.append(trans)
##    entries = list(entries.items())
##    entries.sort()
##    with open("eatT1.txt", 'w', encoding='utf8') as file:
##        for eng, entry in entries:
##            print(ENTRY_SEP, file=file)
##            print("en{}{}{}".format(LANG_ITEM_SEP, eng, LANG_SEP), file=file)
##            other_l = ["{}{}{}".format(lg, LANG_ITEM_SEP, ITEM_SEP.join(trans)) for lg, trans in entry.items()]
##            other_l = (LANG_SEP + "\n").join(other_l)
##            print(other_l, file=file)
##            
##    return entries

##def new_te(write=False):
##    res = []
##    with open("../LingData/Te/tmp.txt", encoding='utf8') as file:
##        for line in file:
##            if ":" in line:
##                root, feats = line.split(':')
##                cls = feats.split('cls=')
##                cls = cls[1].split(',')[0]
##                rc = "{} [cls={}]".format(root, cls)
##                if rc not in res:
##                    res.append(rc)
##    res.sort()
##    if write:
##        with open("../LingData/Te/new_vroots.txt", 'w', encoding='utf8') as file:
##            for r in res:
##                print(r, file=file)
##    return res

##def eng_ti_am(write=True):
##    eng = {}
##    n = 0
##    with open("en_am_v1.txt", encoding='utf8') as file:
##        n += 1
##        if n % 50 == 0:
##            print("Processed {} lines".format(n))
##        for line in file:
##            e, a = line.split(';;')
##            arf, agz = a.split(';')
##            agz = agz.strip()
##            entry = "{}:{}".format(arf, agz)
##            if e in eng:
##                eng[e]['am'].append(entry)
##            else:
##                eng[e] = {'am': [entry]}
##    with open("ti_en_v1.txt", encoding='utf8') as file:
##        for line in file:
##            t, e = line.split(';;')
##            tsplit = t.split(';')
##            tgz, trf = t.split(';')
##            tgz = tgz.strip()
##            tr, tf = trf.split(':')
##            geez = ti_gen(tr, tf, True)
##            geez = geez[-1]
##            entry = "{}:{}".format(trf, geez)
##            for ee in e.strip().split(';'):
##                if ee in eng:
##                    if 'ti' in eng[ee]:
##                        eng[ee]['ti'].append(entry)
##                    else:
##                        eng[ee]['ti'] = [entry]
##                else:
##                    eng[ee] = {'ti': [entry]}
##    if write:
##        with open("eat1.txt", 'w', encoding='utf8') as file:
##            for e, at in eng.items():
##                print(ENTRY_SEP, file=file)
##                print("en;{};;".format(e), file=file)
##                other_l = ["{};{}".format(lg, ';'.join(trans)) for lg, trans in at.items()]
##                other_l = ";;\n".join(other_l)
##                print(other_l, file=file)
###                for lg, trans in at.items():
###                    print("{}:{}".format(lg, trans), file=file)
##    return eng

def proc_en_te_ti():
    results = []
    with open("Tet1.txt", encoding='utf8') as file:
        for line in file:
            te, en, ti = line.split(';;')
            ti = ti.strip()
            tigeez, tirf = ti.split(';')
            tir, tif = tirf.split(":")
            tigeezgem = ti_gen(tir, tif, True)
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

def load(language):
    global AMgen
    global TIgen
    global KSgen
    global geezify
    if language == 'am' and not AMgen:
        AMgen = hm.morpho.get_language('am', segment=False, phon=True).morphology['v']
    elif language == 'ti' and not TIgen:
        TIgen = hm.morpho.get_language('ti', phon=True).morphology['v']
    elif language == 'ks' and not KSgen:
        KSgen = hm.morpho.get_language('gru').morphology['v']
    if not geezify:
        geezify = hm.morpho.geez.geezify

def get_ti_roots():
    roots = []
    with open("../HornMorpho/hm/languages/ti/lex/vb_root.lex") as file:
        for line in file:
            root = line.split()[0]
            if root in roots:
                print("Duplicate: {}".format(root))
            else:
                roots.append(root)
    return roots

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

def ks_gen(root, feats, cls='A', gemination=True):
    load('ks')
    basic = FS("[cls={},sg=m]".format(cls))
    feats = FS(feats)
    basic.update(feats)
#    print("basic {}".format(basic.__repr__()))
    output = KSgen.gen(root, update_feats=basic)
    if output:
        rom = output[0][0]
        geez = geezify(rom, 'gru', gemination=gemination)
        return rom, geez
    else:
        print("Couldn't generate Ks {}, {}".format(root, feats.__repr__()))

def am_gen(root, feats, gemination=False):
    load('am')
    basic = FS(AMTIfeats)
    feats = FS(feats)
    basic.update(feats)
    if 'ob' in feats:
        basic.update(FS("[ob=[-fem]]"))
    output = AMgen.gen(root, update_feats=basic, phon=True)
    if output:
        rom = output[0][0]
        geez = geezify(rom, gemination=gemination)
        return rom, geez
    else:
        print("Couldn't generate Am {}, {}".format(root, feats.__repr__()))

def ti_gen(root, feats, gemination=False):
    load('ti')
    basic = FS(AMTIfeats)
    feats = FS(feats)
    basic.update(feats)
    if 'ob' in feats:
        basic.update(FS("[ob=[-fem]]"))
    output = TIgen.gen(root, 'ti', update_feats=basic, phon=True)
    if output:
        rom = output[0][0]
        geez = geezify(rom, 'ti', gemination=gemination)
        return rom, geez
    else:
        print("Couldn't generate Ti {}, {}".format(root, feats.__repr__()))

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

##def am_verbs(write=True):
##    verbs = []
##    pss = eng_pp()
##    n = 0
##    with open("en_am_v0.txt", encoding='utf8') as file:
##        for line in file:
##            n += 1
##            if n % 100 == 0:
##                print("Processed {} verbs".format(n))
##            en, am = line.strip().split(' || ')
##            root, feats = am.split('_v')
##            forms = am_gen(root, feats, gemination=True)
##            verbs.append((en, root, feats, forms))
##            ps = eng_passive(en, pss)
##            if ps:
##                verbs.append((ps, root, feats, forms))
##    verbs.sort()
##    if write:
##        with open("en_am_v1.txt", 'w', encoding='utf8') as file:
##            for e, r, f, a in verbs:
##                if a:
##                    print("{};;{}:{};{}".format(e, r, f, a[1]), file=file)
##    return verbs

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

##def ti_verbs(write=True):
##    verbs = []
##    n = 0
##    known = 0
##    unknown = 0
##    kverbs = []
##    uverbs = []
##    estems = eng_verbs()
##    epps = eng_pp()
##    with open("ti_en_v0.txt", encoding='utf8') as file:
##        for line in file:
##            n += 1
##            if n % 100 == 0:
##                print("Processed {} lines".format(n))
##            if len(line.split('-')) != 2:
##                print(line)
##            ti, en = line.split("-")
##            ti_ = ti.strip().replace("'", "")
##            anal = ti_anal(ti_)
##            if not en.strip():
##                print("line {}".format(line))
##            en = proc_ti_eng(en.strip(), estems, epps)
##            if not anal:
###                print("No analysis for {}".format(ti_))
##                unknown += 1
##                anal = ti_anal(ti_, guess=True)
##                if not anal:
##                    print("No guess analysis for {}".format(ti_))
##                else:
##                    for a in anal:
##                        uverbs.append("{};{};;{}".format(ti_, a, en))
##            else:
##                known += 1
##                # use only the first analysis
##                anal = anal[0]
##                kverbs.append("{};{};;{}".format(ti_, anal, en))
##    if write:
##        with open("ti_en_v1.txt", 'w', encoding='utf8') as file:
##            for k in kverbs:
##                print(k, file=file)
##            print('###', file=file)
##            for u in uverbs:
##                print(u, file=file)
##    else:
##        return kverbs, uverbs

def eng_passive(defn, pps=None):
    pps = pps or eng_pp()
    dsplit = defn.split()
    if dsplit[0] == "be_v":
        if len(dsplit) > 1 and dsplit[1] in pps:
            ppstem = pps[dsplit[1]]
            return ppstem + "_v[ps]"
    return ''

##def proc_ti_eng(eng, stems=None, pps=None):
##    stems = stems or eng_verbs()
##    pps = pps or eng_pp()
##    result = []
##    eng = eng.split(',')
##    for e in eng:
##        res = ''
##        e = e.strip()
##        e_ = e.split()
##        if not e_:
##            continue
##        e0 = e_[0]
##        if e0 in stems:
##            res0 = e0 + "_v"
##            res = res0 if len(e_) == 1 else ' '.join([res0] + e_[1:])
##            result.append(res)
##            if e0 == 'be' and len(e_) > 1 and e_[1] in pps:
##                stem = pps[e_[1]]
##                result.append(stem + "_v[ps]")
##        else:
##            print("{} not a verb".format(e0))
##    return ';'.join(result)
