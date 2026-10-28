#!/usr/bin/env python3
import random as lrand
import rpy2.robjects as robjects
import argparse
import numpy
import os,sys,math,pickle

def init():
    lrand.seed(1982)
    robjects.r('library(splines)')
    robjects.r('library(stats4)')
    robjects.r('library(survival)')
    robjects.r('library(mvtnorm)')
    robjects.r('library(modeltools)')
    robjects.r('library(coin)')
    robjects.r('library(MASS)')

def get_class_means(class_sl,feats):
    means = {}
    clk = list(class_sl.keys())
    for fk,f in feats.items():
        means[fk] = [numpy.mean((f[class_sl[k][0]:class_sl[k][1]])) for k in clk]
    return clk,means

def save_res(res,filename):
    with open(filename, 'w') as out:
        for k,v in res['cls_means'].items():
            out.write(k+"\t"+str(math.log(max(max(v),1.0),10.0))+"\t")
            if k in res['lda_res_th']:
                for i,vv in enumerate(v):
                    if vv == max(v):
                        out.write(str(res['cls_means_kord'][i])+"\t")
                        break
                out.write(str(res['lda_res'][k]))
            else: out.write("\t")
            out.write( "\t" + (res['wilcox_res'][k] if 'wilcox_res' in res and k in res['wilcox_res'] else "-")+"\n")

def load_data(input_file, nnorm = False):
    with open(input_file, 'rb') as inputf:
        inp = pickle.load(inputf)
    if nnorm: return inp['feats'],inp['cls'],inp['class_sl'],inp['subclass_sl'],inp['class_hierarchy'],inp['norm']
    else: return inp['feats'],inp['cls'],inp['class_sl'],inp['subclass_sl'],inp['class_hierarchy']

def load_res(input_file):
    with open(input_file, 'rb') as inputf:
        inp = pickle.load(inputf)
    return inp['res'],inp['params'],inp['class_sl'],inp['subclass_sl']


def test_kw_r(cls,feats,p,factors):
    robjects.globalenv["y"] = robjects.FloatVector(feats)
    for i,f in enumerate(factors):
        robjects.globalenv['x'+str(i+1)] = robjects.FactorVector(robjects.StrVector(cls[f]))
    fo = "y~x1"
    #for i,f in enumerate(factors[1:]):
    #   if f == "subclass" and len(set(cls[f])) <= len(set(cls["class"])): continue
    #   if len(set(cls[f])) == len(cls[f]): continue
    #   fo += "+x"+str(i+2)
    kw_res = robjects.r('kruskal.test('+fo+',)$p.value')
    return float(tuple(kw_res)[0]) < p, float(tuple(kw_res)[0])

def test_rep_wilcoxon_r(sl,cl_hie,feats,th,multiclass_strat,mul_cor,fn,min_c,comp_only_same_subcl,curv=False):
    comp_all_sub = not comp_only_same_subcl
    tot_ok =  0
    alpha_mtc = th
    all_diff = []
    for pair in [(x,y) for x in cl_hie.keys() for y in cl_hie.keys() if x < y]:
        dir_cmp = "not_set" #
        l_subcl1, l_subcl2 = (len(cl_hie[pair[0]]), len(cl_hie[pair[1]]))
        if mul_cor != 0: alpha_mtc = th*l_subcl1*l_subcl2 if mul_cor == 2 else 1.0-math.pow(1.0-th,l_subcl1*l_subcl2)
        ok = 0
        curv_sign = 0
        first = True
        for i,k1 in enumerate(cl_hie[pair[0]]):
            br = False
            for j,k2 in enumerate(cl_hie[pair[1]]):
                if not comp_all_sub and k1[len(pair[0]):] != k2[len(pair[1]):]:
                    ok += 1
                    continue
                cl1 = feats[sl[k1][0]:sl[k1][1]]
                cl2 = feats[sl[k2][0]:sl[k2][1]]
                med_comp = False
                if len(cl1) < min_c or len(cl2) < min_c:
                    med_comp = True
                sx,sy = numpy.median(cl1),numpy.median(cl2)
                if cl1[0] == cl2[0] and len(set(cl1)) == 1 and  len(set(cl2)) == 1:
                    tres, first = False, False
                elif not med_comp:
                    robjects.globalenv["x"] = robjects.FloatVector(cl1+cl2)
                    robjects.globalenv["y"] = robjects.FactorVector(robjects.StrVector(["a" for a in cl1]+["b" for b in cl2]))
                    pv = float(robjects.r('pvalue(wilcox_test(x~y,data=data.frame(x,y)))')[0])
                    tres = pv < alpha_mtc*2.0
                if first:
                    first = False
                    if not curv and ( med_comp or tres ):
                        dir_cmp = sx < sy
                        #if sx == sy: br = True
                    elif curv:
                        dir_cmp = None
                        if med_comp or tres:
                            curv_sign += 1
                            dir_cmp = sx < sy
                    else: br = True
                elif not curv and med_comp:
                    if ((sx < sy) != dir_cmp or sx == sy): br = True
                elif curv:
                    if tres and dir_cmp == None:
                        curv_sign += 1
                        dir_cmp = sx < sy
                    if tres and dir_cmp != (sx < sy):
                        br = True
                        curv_sign = -1
                elif not tres or (sx < sy) != dir_cmp or sx == sy: br = True
                if br: break
                ok += 1
            if br: break
        if curv: diff = curv_sign > 0
        else: diff = (ok == len(cl_hie[pair[1]])*len(cl_hie[pair[0]])) # or (not comp_all_sub and dir_cmp != "not_set")
        if diff: tot_ok += 1
        if not diff and multiclass_strat: return False
        if diff and not multiclass_strat: all_diff.append(pair)
    if not multiclass_strat:
        tot_k = len(list(cl_hie.keys()))
        for k in cl_hie.keys():
            nk = 0
            for a in all_diff:
                if k in a: nk += 1
            if nk == tot_k-1: return True
        return False
    return True



def contast_within_classes_or_few_per_class(feats,inds,min_cl,ncl):
    ff = list(zip(*[v for n,v in feats.items() if n != 'class']))
    cols = [ff[i] for i in inds]
    cls = [feats['class'][i] for i in inds]
    if len(set(cls)) < ncl:
        return True
    for c in set(cls):
        if cls.count(c) < min_cl:
            return True
        cols_cl = [x for i,x in enumerate(cols) if cls[i] == c]
        for i,col in enumerate(zip(*cols_cl)):
            if (len(set(col)) <= min_cl and min_cl > 1) or (min_cl == 1 and len(set(col)) <= 1):
                return True
    return False

# Local addition, not upstream: True if R's MASS::lda would stop on this
# subsample, i.e. it holds a single class, or some feature is constant within
# every class (zero pooled within-class variance).
def lda_would_fail(feats,inds):
    cls = [feats['class'][i] for i in inds]
    if len(set(cls)) < 2:
        return True
    for k,v in feats.items():
        if k in ('class','subclass','subject'):
            continue
        if all(len(set(v[i] for i,ci in zip(inds,cls) if ci == c)) <= 1 for c in set(cls)):
            return True
    return False

def test_lda_r(cls,feats,cl_sl,boots,fract_sample,lda_th,tol_min,nlogs):
    fk = list(feats.keys())
    means = dict([(k,[]) for k in feats.keys()])
    feats['class'] = list(cls['class'])
    clss = list(set(feats['class']))

    for uu,k in enumerate(fk):
        if k == 'class':
            continue

        ff = [(feats['class'][i],v) for i,v in enumerate(feats[k])]

        for c in clss:
            if len(set([float(v[1]) for v in ff if v[0] == c])) > max(float(feats['class'].count(c))*0.5,4):
                continue

            for i,v in enumerate(feats[k]):
                if feats['class'][i] == c:
                    feats[k][i] = math.fabs(feats[k][i] + lrand.normalvariate(0.0,max(feats[k][i]*0.05,0.01)))

    rdict = {}

    for a,b in feats.items():
        if a == 'class' or a == 'subclass' or a == 'subject':
            rdict[a] = robjects.StrVector(b)
        else:
            rdict[a] = robjects.FloatVector(b)

    robjects.globalenv["d"] = robjects.DataFrame(rdict)
    lfk = len(feats[fk[0]])
    rfk = int(float(len(feats[fk[0]]))*fract_sample)
    f = "class ~ "+fk[0]

    for k in fk[1:]:
        f += " + " + k.strip()

    ncl = len(set(cls['class']))
    min_cl = int(float(min([cls['class'].count(c) for c in set(cls['class'])]))*fract_sample*fract_sample*0.5)
    min_cl = max(min_cl,1)
    pairs = [(a,b) for a in set(cls['class']) for b in set(cls['class']) if a > b]

    for k in fk:
        for i in range(boots):
            means[k].append([])

    # Local change, not upstream: upstream silently keeps the last draw when
    # none of the 1000 passes the per-class check, and R's lda then stops if
    # that draw is degenerate. Here that case is counted and warned about, and
    # a degenerate draw is redrawn until lda can fit it. A last draw that lda
    # can fit is kept, so runs that worked upstream are unchanged.
    n_unchecked = n_redrawn = 0
    for i in range(boots):
        for rtmp in range(1000):
            rand_s = [lrand.randint(0,lfk-1) for v in range(rfk)]
            if not contast_within_classes_or_few_per_class(feats,rand_s,min_cl,ncl):
                break
        else:
            n_unchecked += 1
            if lda_would_fail(feats,rand_s):
                n_redrawn += 1
                for rtmp in range(1000):
                    rand_s = [lrand.randint(0,lfk-1) for v in range(rfk)]
                    if not lda_would_fail(feats,rand_s):
                        break
                else:
                    sys.exit("LEfSe LDA: no subsample of "+str(rfk)+" samples that lda can fit was found in 1000 draws; this design is too small for LEfSe's LDA step.")

        rand_s = [r+1 for r in rand_s]
        means[k][i] = []

        for p in pairs:
            robjects.globalenv["rand_s"] = robjects.IntVector(rand_s)
            robjects.globalenv["sub_d"] = robjects.r('d[rand_s,]')
            z = robjects.r('z <- suppressWarnings(lda(as.formula('+f+'),data=sub_d,tol='+str(tol_min)+'))')
            robjects.r('w <- z$scaling[,1]')
            robjects.r('w.unit <- w/sqrt(sum(w^2))')
            robjects.r('ss <- sub_d[,-match("class",colnames(sub_d))]')

            if 'subclass' in feats:
                robjects.r('ss <- ss[,-match("subclass",colnames(ss))]')

            if 'subject' in feats:
                robjects.r('ss <- ss[,-match("subject",colnames(ss))]')

            robjects.r('xy.matrix <- as.matrix(ss)')
            robjects.r('LD <- xy.matrix%*%w.unit')
            robjects.r('effect.size <- abs(mean(LD[sub_d[,"class"]=="'+p[0]+'"]) - mean(LD[sub_d[,"class"]=="'+p[1]+'"]))')
            scal = robjects.r('wfinal <- w.unit * effect.size')
            rres = robjects.r('mm <- z$means')
            rowns = list(rres.rownames)
            lenc = len(list(rres.colnames))
            coeff = [abs(float(v)) if not math.isnan(float(v)) else 0.0 for v in scal]
            res = dict([(pp,[float(ff) for ff in rres.rx(pp,True)] if pp in rowns else [0.0]*lenc ) for pp in [p[0],p[1]]])

            for j,k in enumerate(fk):
                gm = abs(res[p[0]][j] - res[p[1]][j])
                means[k][i].append((gm+coeff[j])*0.5)

    if n_unchecked:
        sys.stderr.write("WARNING: LEfSe LDA: in "+str(n_unchecked)+" of "+str(boots)+" bootstrap iterations no subsample of "+str(rfk)+" samples passed LEfSe's per-class check in 1000 draws, so LDA ran on an unchecked subsample ("+str(n_redrawn)+" redrawn because lda could not fit the first one). LDA scores rest on weaker subsamples than LEfSe intends.\n")

    res = {}

    for k in fk:
        m = max([numpy.mean([means[k][kk][p] for kk in range(boots)]) for p in range(len(pairs))])
        res[k] = math.copysign(1.0,m)*math.log(1.0+math.fabs(m),10)

    return res,dict([(k,x) for k,x in res.items() if math.fabs(x) > lda_th])


def test_svm(cls,feats,cl_sl,boots,fract_sample,lda_th,tol_min,nsvm):
    return None
"""
    fk = feats.keys()
    clss = list(set(cls['class']))
    y = [clss.index(c)*2-1 for c in list(cls['class'])]
    xx = [feats[f] for f in fk]
    if nsvm:
        maxs = [max(v) for v in xx]
        mins = [min(v) for v in xx]
        x = [ dict([(i+1,(v-mins[i])/(maxs[i]-mins[i])) for i,v in enumerate(f)]) for f in zip(*xx)]
    else: x = [ dict([(i+1,v) for i,v in enumerate(f)]) for f in zip(*xx)]

    lfk = len(feats[fk[0]])
    rfk = int(float(len(feats[fk[0]]))*fract_sample)
    mm = []

    best_c = svmutil.svm_ms(y, x, [pow(2.0,i) for i in range(-5,10)],'-t 0 -q')
    for i in range(boots):
        rand_s = [lrand.randint(0,lfk-1) for v in range(rfk)]
        r = svmutil.svm_w([y[yi] for yi in rand_s], [x[xi] for xi in rand_s], best_c,'-t 0 -q')
        mm.append(r[:len(fk)])
    m = [numpy.mean(v) for v in zip(*mm)]
    res = dict([(v,m[i]) for i,v in enumerate(fk)])
    return res,dict([(k,x) for k,x in res.items() if math.fabs(x) > lda_th])
"""

def read_params(args):
    parser = argparse.ArgumentParser(description='LEfSe 1.1.01')
    parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="the input file")
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str,
                help="the output file containing the data for the visualization module")
    parser.add_argument('-o',dest="out_text_file", metavar='str', type=str, default="",
                help="set the file for exporting the result (only concise textual form)")
    parser.add_argument('-a',dest="anova_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Anova test (default 0.05)")
    parser.add_argument('-w',dest="wilcoxon_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Wilcoxon test (default 0.05)")
    parser.add_argument('-l',dest="lda_abs_th", metavar='float', type=float, default=2.0,
                help="set the threshold on the absolute value of the logarithmic LDA score (default 2.0)")
    parser.add_argument('--nlogs',dest="nlogs", metavar='int', type=int, default=3,
        help="max log ingluence of LDA coeff")
    parser.add_argument('--verbose',dest="verbose", metavar='int', choices=[0,1], type=int, default=0,
        help="verbose execution (default 0)")
    parser.add_argument('--wilc',dest="wilc", metavar='int', choices=[0,1], type=int, default=1,
        help="wheter to perform the Wicoxon step (default 1)")
    parser.add_argument('-r',dest="rank_tec", metavar='str', choices=['lda','svm'], type=str, default='lda',
        help="select LDA or SVM for effect size (default LDA)")
    parser.add_argument('--svm_norm',dest="svm_norm", metavar='int', choices=[0,1], type=int, default=1,
        help="whether to normalize the data in [0,1] for SVM feature waiting (default 1 strongly suggested)")
    parser.add_argument('-b',dest="n_boots", metavar='int', type=int, default=30,
                help="set the number of bootstrap iteration for LDA (default 30)")
    parser.add_argument('-e',dest="only_same_subcl", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test only among the subclasses with the same name (default 0)")
    parser.add_argument('-c',dest="curv", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test ing the Curtis's approach [BETA VERSION] (default 0)")
    parser.add_argument('-f',dest="f_boots", metavar='float', type=float, default=0.67,
                help="set the subsampling fraction value for each bootstrap iteration (default 0.66666)")
    parser.add_argument('-s',dest="strict", choices=[0,1,2], type=int, default=0,
                help="set the multiple testing correction options. 0 no correction (more strict, default), 1 correction for independent comparisons, 2 correction for dependent comparison")
#       parser.add_argument('-m',dest="m_boots", type=int, default=5,
#               help="minimum cardinality of classes in each bootstrapping iteration (default 5)")
    parser.add_argument('--min_c',dest="min_c", metavar='int', type=int, default=10,
                help="minimum number of samples per subclass for performing wilcoxon test (default 10)")
    parser.add_argument('-t',dest="title", metavar='str', type=str, default="",
                help="set the title of the analysis (default input file without extension)")
    parser.add_argument('-y',dest="multiclass_strat", choices=[0,1], type=int, default=0,
                help="(for multiclass tasks) set whether the test is performed in a one-against-one ( 1 - more strict!) or in a one-against-all setting ( 0 - less strict) (default 0)")
    args = parser.parse_args()

    params = vars(args)
    if params['title'] == "":
        params['title'] = params['input_file'].split("/")[-1].split('.')[0]

    return params


def lefse_run():
    init()
    params = read_params(sys.argv)
    feats,cls,class_sl,subclass_sl,class_hierarchy = load_data(params['input_file'])
    kord,cls_means = get_class_means(class_sl,feats)
    wilcoxon_res = {}
    kw_n_ok = 0
    nf = 0
    for feat_name,feat_values in list(feats.items()):
        if params['verbose']:
            print("Testing feature",str(nf),": ",feat_name)
            nf += 1
        kw_ok,pv = test_kw_r(cls,feat_values,params['anova_alpha'],sorted(cls.keys()))
        if not kw_ok:
            if params['verbose']: print("\tkw ko")
            del feats[feat_name]
            wilcoxon_res[feat_name] = "-"
            continue
        if params['verbose']: print("\tkw ok\t")

        if not params['wilc']: continue
        kw_n_ok += 1
        res_wilcoxon_rep = test_rep_wilcoxon_r(subclass_sl,class_hierarchy,feat_values,params['wilcoxon_alpha'],params['multiclass_strat'],params['strict'],feat_name,params['min_c'],params['only_same_subcl'],params['curv'])
        wilcoxon_res[feat_name] = str(pv) if res_wilcoxon_rep else "-"
        if not res_wilcoxon_rep:
            if params['verbose']: print("wilc ko")
            del feats[feat_name]
        elif params['verbose']: print("wilc ok\t")

    if len(feats) > 0:
        print("Number of significantly discriminative features:", len(feats), "(", kw_n_ok, ") before internal wilcoxon")
        if params['lda_abs_th'] < 0.0:
            lda_res,lda_res_th = dict([(k,0.0) for k,v in feats.items()]), dict([(k,v) for k,v in feats.items()])
        else:
            if params['rank_tec'] == 'lda': lda_res,lda_res_th = test_lda_r(cls,feats,class_sl,params['n_boots'],params['f_boots'],params['lda_abs_th'],0.0000000001,params['nlogs'])
            elif params['rank_tec'] == 'svm': lda_res,lda_res_th = test_svm(cls,feats,class_sl,params['n_boots'],params['f_boots'],params['lda_abs_th'],0.0,params['svm_norm'])
            else: lda_res,lda_res_th = dict([(k,0.0) for k,v in feats.items()]), dict([(k,v) for k,v in feats.items()])
    else:
        print("Number of significantly discriminative features:", len(feats), "(", kw_n_ok, ") before internal wilcoxon")
        print("No features with significant differences between the two classes")
        lda_res,lda_res_th = {},{}
    outres = {}
    outres['lda_res_th'] = lda_res_th
    outres['lda_res'] = lda_res
    outres['cls_means'] = cls_means
    outres['cls_means_kord'] = kord
    outres['wilcox_res'] = wilcoxon_res
    print("Number of discriminative features with abs LDA score >",params['lda_abs_th'],":",len(lda_res_th))
    save_res(outres,params["output_file"])


if __name__ == '__main__':
    lefse_run()