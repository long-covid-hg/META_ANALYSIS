#!/usr/bin/env python3
import argparse
import json
import gzip
from collections import namedtuple
import sys
import math
from scipy.stats import chi2
import scipy.stats
import numpy
from typing import Dict, Tuple, List
import subprocess
from collections import deque
import re

chrord = { "chr"+str(chr):int(chr) for chr in list(range(1,23))}
chrord["X"] = 23
chrord["Y"] = 24
chrord["MT"] = 25
chrord["chrX"] = 23
chrord["chrY"] = 24
chrord["chrMT"] = 25
chrord.update({str(chr):int(chr) for chr in list(range(1,25)) } )
chrord.update({int(chr):int(chr) for chr in list(range(1,25)) } )

re_allele = re.compile('^[ATCG]+$', re.IGNORECASE)


def het_test( effs_sizes: List[float], weights: List[float], effs_size_meta: float) -> float:
    '''
        Computes Cochran's Q test for heterogeneity
        input:
            effs_sizes: original effect sizes
            weights: weights
            effs_size_meta: effect size from meta-analysis
        output:
            p-value
    '''
    k=len(effs_sizes)

    effs_sizes_array=numpy.array(effs_sizes)
    weights_array=numpy.array(weights)
    eff_dev=weights_array*((effs_sizes_array-effs_size_meta)**2)
    sum_eff_dev=numpy.sum(eff_dev)

    return scipy.stats.distributions.chi2.sf(sum_eff_dev, k-1)

def n_meta( studies : List[Tuple['Study','VariantData']] ) -> Tuple:
    '''
        Computes sample size weighted meta-analysis for variants in studies
        input:
            studies: studies and data in tuples
        output:
            tuple with results from meta-analysis or None
    '''
    weights = []
    effs_size_org = []

    effs_size = []
    tot_size =0
    sum_betas=0
    sum_weights=0

    sum_af_alt=0
    sum_info=0
    sum_Nsamples=0

    for s in studies:
        study = s[0]
        dat = s[1]

        if dat.Nsamples is None or dat.Nsamples==0:
            print("Number of samples was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            return None

        effs_size.append( math.sqrt(study.effective_size) * numpy.sign(dat.beta) * dat.z_score)
        sum_weights+= math.sqrt(study.effective_size)
        sum_betas+=math.sqrt(study.effective_size) * dat.beta
        tot_size+=study.effective_size
        weights.append(math.sqrt(study.effective_size))
        effs_size_org.append(dat.beta)

        sum_af_alt+=(dat.af_alt*dat.Nsamples)
        sum_info+=(dat.info*dat.Nsamples)
        sum_Nsamples+=(dat.Nsamples)

    beta_meta=sum_betas/sum_weights

    af_alt_meta=sum_af_alt/sum_Nsamples
    info_meta=sum_info/sum_Nsamples


    #TODO se
    return ( beta_meta, None, max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf( abs( sum( effs_size ) ) / math.sqrt(tot_size) )), af_alt_meta, info_meta, effs_size_org, weights) if len(effs_size)==len(studies) else None


def inv_var_meta( studies : List[Tuple['Study','VariantData']] ) -> Tuple:
    '''
        Computes inverse-variance weighted meta-analysis for variants in studies
        input:
            studies: studies and data in tuples
        output:
            tuple with results from meta-analysis or None
    '''
    weights = []
    effs_size_org = []

    effs_inv_var = []
    sum_inv_var=0
    sum_af_alt=0
    sum_info=0
    sum_Nsamples=0

    for s in studies:
        study = s[0]
        dat = s[1]
        if dat.se is None or dat.se==0:
            print("Standard error was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            return None
        # NEED TO USE n_cases AND n_controls IF Nsamples NOT PROVIDED...
        if dat.Nsamples is None or dat.Nsamples==0:
            print("Number of samples was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            return None

        var = (dat.se * dat.se)

        inv_var =  (1/var)
        sum_inv_var+=inv_var
        effs_inv_var.append( inv_var *  dat.beta )

        weights.append(inv_var)
        effs_size_org.append(dat.beta)

        sum_af_alt+=(dat.af_alt*dat.Nsamples)
        sum_info+=(dat.info*dat.Nsamples)
        sum_Nsamples+=(dat.Nsamples)

    beta_meta=sum(effs_inv_var)/ sum_inv_var

    af_alt_meta=sum_af_alt/sum_Nsamples
    info_meta=sum_info/sum_Nsamples

    return (beta_meta, math.sqrt(1/sum_inv_var), max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf(abs(sum(effs_inv_var) / math.sqrt(sum_inv_var) ))), af_alt_meta, info_meta, effs_size_org, weights)


def variance_weight_meta( studies : List[Tuple['Study','VariantData']] ) -> Tuple:
    '''
        Computes variance weighted meta-analysis for variants in studies
        input:
            studies: studies and data in tuples
        output:
            tuple with results from meta-analysis or None
    '''
    weights = []
    effs_size_org = []

    effs_se = []
    tot_se = 0
    sum_weights=0
    sum_betas=0
    for s in studies:
        study = s[0]
        dat = s[1]

        if dat.se is None or dat.se==0:
            print("Standard error was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            return None
        weight =  (1/dat.se) * dat.z_score
        sum_weights+=weight
        sum_betas+= weight * dat.beta
        effs_se.append( weight * numpy.sign(dat.beta)  )
        tot_se+=1/ (dat.se * dat.se)

        weights.append(weight)
        effs_size_org.append(dat.beta)

    beta_meta=sum_betas / sum_weights

    #TODO SE
    return (beta_meta, None, max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf( abs( sum( effs_se ) ) /  math.sqrt(tot_se))), effs_size_org, weights)


SUPPORTED_METHODS = {"n":n_meta,"inv_var":inv_var_meta,"variance":variance_weight_meta}


def check_eff_field(field):
    if field.lower() in ["beta","or"]:
        return field.lower()
    else:
        raise Exception("effect_type must be beta or OR")

flip = {"A":"T","C":"G","T":"A","G":"C"}

def flip_strand( allele):
    return "".join([ flip[a] for a in allele])

def is_symmetric(a1, a2):
    return (a1=="A" and a2=="T") or (a1=="T" and a2=="A") or (a1=="C" and a2=="G") or (a1=="G" and a2=="C")


class VariantData:

    def __init__(self, chr, pos, ref, alt, af_alt, info, beta, pval, se=None, Nsamples=None, extra_cols=[]):
        self.chr = chr
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af_alt = af_alt
        self.info = info
        self.beta = beta
        self.pval = pval
        try:
            self.se = float(se) if se is not None else None
            self.Nsamples = int(float(Nsamples)) if Nsamples is not None else None
        except ValueError:
            self.se = None
            self.Nsamples = None

        self.extra_cols = extra_cols

    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __lt__(self, other):

        return (  (self.chr==other.chr and self.pos<other.pos)
                  or (self.chr < other.chr)
               )

    def is_equal(self, other:'VariantData') -> bool:
        """
            Checks if this VariantData is the same variant (possibly different strand or ordering of alleles)
            returns: true if the same false if not
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if self.ref== other.ref and self.alt == other.alt :
                return True

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                if self.ref == other.ref and self.alt == other.alt:
                    return True
                elif self.ref == other.alt and self.alt == other.ref:
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                return True

        return False

    def equalize_to(self, other:'VariantData') -> bool:
        """
            Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles)
            If it is, changes this variant's alleles and beta accordingly
            returns: true if the same (flips effect direction and ref/alt alleles if necessary) or false if not the same variant
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if self.ref== other.ref and self.alt == other.alt :
                    return True

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                if self.ref == other.ref and self.alt == other.alt:
                    return True
                elif self.ref == other.alt and self.alt == other.ref:
                    self.beta = -1 * self.beta if self.beta is not None else None
                    self.af_alt = 1 - self.af_alt if self.af_alt is not None else None
                    t = self.alt
                    self.alt = self.ref
                    self.ref = t
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af_alt = 1 - self.af_alt if self.af_alt is not None else None
                t = self.alt
                self.alt = self.ref
                self.ref = t
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                self.ref = flip_strand(self.ref)
                self.alt = flip_strand(self.alt)
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af_alt = 1 - self.af_alt if self.af_alt is not None else None
                self.ref =flip_strand(self.alt)
                self.alt = flip_strand(self.ref)
                return True

        return False

    @property
    def z_score(self):
        '''
            Lazy compute unsigned z-score
        '''
        if self.z_scr is None:
            self.z_scr = math.sqrt(chi2.isf(self.pval, df=1))
        return self.z_scr

    def __str__(self):
        return "chr:{} pos:{} ref:{} alt:{} beta:{} pval:{} se:{} ".format(self.chr, self.pos, self.ref, self.alt, self.beta, self.pval, self.se)


class Study:
    REQUIRED_DATA_FIELDS = {"chr":str,"pos":str,"ref":str,"alt":str,"af_alt":str,"info":str,"effect":str,"pval":str}

    REQUIRED_CONF = {"name":str,"file":str, "n_cases": int, "n_controls":int,
    "chr":str,"pos":str,"ref":str,"alt":str,"af_alt":str,"info":str,"effect":str,
    "effect_type":check_eff_field,
    "pval":str}

    OPTIONAL_FIELDS = {"se":str,"Nsamples":str}

    def __init__(self, conf, chrom=None, dont_allow_space=False):
        '''
        chrom: a chromosome to limit to or None if all chromosomes
        dont_allow_space: boolean, don't treat space as field delimiter (only tab)
        '''
        self.conf =conf
        self.chrom = chrord[chrom] if chrom is not None else None
        self.dont_allow_space = dont_allow_space
        self.future = deque()
        self.eff_size= None
        self.eff_size_logistic= None
        self.z_scr = None
        self.prev_var = None
        for v in Study.REQUIRED_CONF:
            if v not in self.conf:
                raise Exception("Meta configuration for study must contain required elements: "
                    + ",".join(Study.REQUIRED_CONF.keys() ) + ". Offending configuration: " + str(self.conf))

            try:
                self.conf[v] = Study.REQUIRED_CONF[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + str(v) +
                    " in configuration: " + str(self.conf) + ". ERR:" + str(e))

        for v in Study.OPTIONAL_FIELDS:
            if v not in self.conf:
                continue
            try:
                self.conf[v] = Study.OPTIONAL_FIELDS[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + v +
                    " in configuration: " + str(self.conf) + ". ERR:" + str(e))

        self.conf["fpoint"] = gzip.open(conf["file"],'rt')
        if self.dont_allow_space:
            header = conf["fpoint"].readline().rstrip().split('\t')
        else:
            header = conf["fpoint"].readline().rstrip().split()

        for k in Study.REQUIRED_DATA_FIELDS.keys():
            if self.conf[k] not in header:
                raise Exception("Required headers not in data in study " + self.conf["name"] + ". Missing:" + ",".join([ self.conf[k] for k in Study.REQUIRED_DATA_FIELDS.keys() if self.conf[k] not in header])  )
        self.conf["h_idx"] = { k:header.index( self.conf[k] ) for k in Study.REQUIRED_DATA_FIELDS.keys() }

        for f in Study.OPTIONAL_FIELDS.keys():
            if f in self.conf:
                 if self.conf[f] not in header:
                     raise Exception("Configured column " + self.conf[f] + " not found in the study results " + self.conf["name"])
                 self.conf["h_idx"][f] = header.index(self.conf[f])

        if "extra_cols" in self.conf:
            for c in self.conf["extra_cols"]:
                if c not in header:
                    raise Exception("Configured column " + self.conf[c] + " not found in the study results " + self.conf["name"])
                self.conf["h_idx"][c] = header.index(c)
        else:
             self.conf["extra_cols"] = []



    @property
    def n_cases(self):
        return self.conf["n_cases"]

    @property
    def n_controls(self):
        return self.conf["n_controls"]

    @property
    def n_samples(self):
        return self.conf["n_cases"]+self.conf["n_controls"]


    @property
    def effective_size(self):
        if self.eff_size is None:
            self.eff_size = ( (4 * self.n_cases *  self.n_controls  ) / ( self.n_cases+  self.n_controls ))
        return self.eff_size
    @property
    def effective_size_logistic(self):
        if self.eff_size_logistic is None:
            self.eff_size_logistic = self.n_cases * (1 - self.n_cases / float(self.n_cases + self.n_controls))
        return self.eff_size_logistic
    @property
    def name(self):
        return self.conf["name"]

    def has_std_err(self):
        return "se" in self.conf

    def has_n_samples(self):
        return "Nsamples" in self.conf

    def get_next_data(self, just_one: bool = False) -> List[VariantData]:
        """
            Returns a list of variants. List containts >1 elements if they are on the same position and just_one == False.

            input:
                just_one: always returns only the next variant in order and not all next with the same position
            output:
                list of next variants
        """

        vars = list()

        if len(self.future)>0:
            ## only return variants with same position so that possible next variant position stored stays
            f = [ (i,v) for i,v in enumerate(self.future) if i==0 or (v.chr==self.future[0].chr and v.pos==self.future[0].pos) ]
            for i,v in reversed(f):
                del self.future[i]
            vars.extend([ v for i,v in f ])
            if len(self.future)>0:
                return vars

        while True:
            chr = None
            ref = None
            alt = None
            l = None
            ## loop ignoring  alternate contigs and non-ATCG alleles for now.
            while chr is None or chr not in chrord or (self.chrom is not None and chr != self.chrom) or re_allele.match(ref) is None or re_allele.match(alt) is None:
                l = self.conf["fpoint"].readline()
                if l=="":
                    return None if len(vars) == 0 else vars

                if self.dont_allow_space:
                    l = l.rstrip().split('\t')
                else:
                    l = l.rstrip().split()
                chr = l[self.conf["h_idx"]["chr"]]
                chr = chrord[chr] if chr in chrord else None
                ref = l[self.conf["h_idx"]["ref"]]
                alt = l[self.conf["h_idx"]["alt"]]

            pos = l[self.conf["h_idx"]["pos"]]
            af_alt = l[self.conf["h_idx"]["af_alt"]]
            info = l[self.conf["h_idx"]["info"]]
            eff = l[self.conf["h_idx"]["effect"]]
            pval = l[self.conf["h_idx"]["pval"]]

            pos = int(float(pos))

            se = l[self.conf["h_idx"]["se"]] if "se" in self.conf["h_idx"] else None
            Nsamples = l[self.conf["h_idx"]["Nsamples"]] if "Nsamples" in self.conf["h_idx"] else None

            effect_type = self.conf["effect_type"]
            try:
                af_alt = float(af_alt)
                info = float(info)
                pval = float(pval)
                eff = float(eff)
            except Exception as e:
                af_alt = None
                info = None
                pval = None
                eff = None

            if( effect_type=="or" and eff):
                eff = math.log(eff)

            extracols = [ l[self.conf["h_idx"][c]] for c in self.conf["extra_cols"] ]

            v = VariantData(chr,pos,ref,alt,af_alt,info,eff,pval,se,Nsamples,extracols)

            if self.prev_var is not None and v < self.prev_var:
                raise Exception("Disorder in study " + self.conf['name'] + " in file " + self.conf['file'] + ". Sort all summary statistic files by chromosome and then position and rerun.\nOffending line: " + "\t".join(l))
            self.prev_var = v
            if len(vars)==0 or ( vars[0].chr == v.chr and vars[0].pos == v.pos  ):
                added=False
                for v_ in vars:
                    if v == v_:
                        print('ALREADY ADDED FOR STUDY ' + self.name + ': ' + str(v), file=sys.stderr)
                        added=True
                if not added:
                    vars.append(v )
                if just_one:
                    break
            else:
                self.future.append(v )
                break

        return vars

    @property
    def extra_cols(self):
        return self.conf["extra_cols"]

    def get_match(self, dat: VariantData) -> VariantData:
        """
            Reads current study until variant in 'dat' is reached or overtaken in chr pos orded.
            IF matching variant found (can flip alleles) the matching VariantData(effect flipped if alleles flipped) is returned.
            input:
                dat: the variant to look for
            output: matching VariantData in this study or None if no match.
        """

        otherdats = self.get_next_data( )

        if otherdats is None or len(otherdats)==0:
            return None

        while otherdats is not None and (otherdats[0].chr<dat.chr or (otherdats[0].chr==dat.chr and otherdats[0].pos<dat.pos)):
            otherdats = self.get_next_data()

        if otherdats is None:
            return None

        if otherdats[0].chr > dat.chr or otherdats[0].pos> dat.pos:
            self.put_back(otherdats)
            return None

        for i,v in enumerate(otherdats):
            if v.equalize_to(dat):
                del otherdats[i]
                self.put_back(otherdats)
                return v

        ## no match but stayed in the same pos. add variants back to future queue
        self.put_back(otherdats)
        return None


    def put_back(self, variantlist: List[VariantData]):
        '''
        Put list of variants back to wait for matching

        input:
            variantlist: list of VariantData objects
        output:
            p-value
        '''

        self.future.extendleft(variantlist)


def get_studies(conf:str, chrom, dont_allow_space) -> List[Study]:
    """
        Reads json configuration and returns studies in the meta
    """

    studies_conf = json.load(open(conf,'r'))
    std_list = studies_conf["meta"]

    return [ Study(s, chrom, dont_allow_space) for s in studies_conf["meta"]]

def do_meta(study_list: List[ Tuple[Study, VariantData]], methods: List[str], is_het_test) -> List[Tuple] :
    '''
        Computes meta-analysis between all studies and data given in the std_list
        input:
            study_list: studies and data in tuples
            methods: list of methods to calculate
            is_het_test: boolean, do heterogeneity test
        output:
            list of tuples (effect_size, standard error, p-value, n_cases, n_controls, n_eff, (, het test p-value)) for each method in the same order as methods were given
    '''
    met = [ SUPPORTED_METHODS[m](study_list) for m in methods ]

    meta_res = []
    n_cases = sum([tuple[0].n_cases for tuple in study_list])
    n_controls = sum([tuple[0].n_controls for tuple in study_list])
    # ADD IN WEIGHTED AF_ALT AND INFO HERE
    n_eff = sum([tuple[0].effective_size_logistic for tuple in study_list])
    for m in met:
        if m is not None:
            if is_het_test:
                meta_res.append((m[0], m[1], m[2], n_cases, n_controls, n_eff, m[3], m[4], het_test(m[5], m[6], m[0])))
            else:
                meta_res.append((m[0], m[1], m[2], n_cases, n_controls, n_eff, m[3], m[4]))
        else:
            meta_res.append(None)

    return meta_res

def format_num(num, precision=5):
    return numpy.format_float_scientific(num, precision=precision) if num is not None else "NA"

def get_next_variant( studies : List[Study]) -> List[VariantData]:
    '''
        get variant data for all studies
        The variant data is the first in chromosomal order across studies (ties broken by alphabetic order of ref)
        input:
            study_list: studies and data in tuples
        output:
            List of VariantData objects. The list is in the same order as the input studies. If smallest variant is
            not found in a study, that position in the list will be Null
    '''

    dats = []
    first = None
    for s in studies:
        d = s.get_next_data()
        dats.append(d)

        if d is not None:
            for v in d:
                if first is None or v < first:
                    first = v

    res = []
    for i,s in enumerate(studies):

        if dats[i] is None:
            res.append(None)
            continue

        # Flag tracks that only one variant (best match) is returned per study, if multiple variants equal to first
        added=False
        for j,v in reversed(list(enumerate(dats[i]))):
            if v == first:
                res.append(v)
                added=True
                del dats[i][j]
                s.put_back(dats[i])
                break
            if not v.is_equal(first):
                s.put_back([v])
                del dats[i][j]
        if not added:
            for j,v in reversed(list(enumerate(dats[i]))):
                if v.equalize_to(first):
                    res.append(v)
                    added=True
                    del dats[i][j]
                    break
            s.put_back(dats[i])
        if not added:
            res.append(None)

    return res



def run():
    '''
        First parameter should be a path to a json configuration file with these elements:
            "name":"FINNGEN",
            "file":"/Users/mitja/projects/finngen/META_ANALYSIS/I9_AF.gz",
            "n_cases": 6570 ,
            "n_controls": 48378,
            "chr":"CHR",
            "pos":"POS",
            "ref":"Allele1",
            "alt":"Allele2",
            "af_alt":"af_alt",
            "info":"info"
            "effect":"BETA",
            "effect_type":"beta",
            "pval":"p.value"
            "se":"SE" <- this parameter is optional. If given for compared studies additional p-value will be added using this as a weight for z-score.
            "Nsamples":"Nsamples" <- this parameter is also optional, as it can be inferred by other parameters.
        Second parameter should be a path to (empty/not existing) directory where the data should be stored
    '''

    parser = argparse.ArgumentParser(description="Run x-way meta-analysis")
    parser.add_argument('config_file', action='store', type=str, help='Configuration file ')
    parser.add_argument('path_to_res', action='store', type=str, help='Result file')

    parser.add_argument('methods', action='store', type=str, help='List of meta-analysis methods to compute separated by commas.'
            + 'Allowed values [n,inv_var,variance]', default="inv_var")

    parser.add_argument('--not_quiet', action='store_false', dest='quiet', help='Print matching variants to stdout')
    parser.set_defaults(quiet=True)

    parser.add_argument('--leave_one_out', action='store_true', help='Do leave-one-out meta-analysis')
    parser.set_defaults(leave_one_out=False)

    parser.add_argument('--leave_most_sig_out', action='store_true', help='Report meta stats with largest study left out')
    parser.set_defaults(leave_most_sig_out=False)

    parser.add_argument('--is_het_test', action='store_true', help='Do heterogeneity tests based on Cochrans Q and output het_p')
    parser.set_defaults(het_test=False)

    parser.add_argument('--pairwise_with_first', action='store_true', help='Do pairwise meta-analysis with the first given study')
    parser.add_argument('--dont_allow_space', action='store_true', help='Do not allow space as field delimiter')

    parser.add_argument('--chrom', action='store', type=str, help='Restrict to given chromosome')

    args = parser.parse_args()

    studs = get_studies(args.config_file, args.chrom, args.dont_allow_space)

    methods = []

    for m in args.methods.split(","):
        if m not in SUPPORTED_METHODS:
            raise Exception("Unsupported meta method" + m + " given. Supported values" + ",".join(SUPPORTED_METHODS))
        methods.append(m)

    if "inv_var" in methods or "variance" in methods:
        for s in studs:
            if not s.has_std_err():
                raise Exception("Variance based method requested but not all studies have se column specified.")
            # THIS SHOULDN'T BE NECESSARY - NEED TO COMBINE n_cases AND n_controls IF Nsamples NOT DEFINED
            if not s.has_n_samples():
                raise Exception("Variance based method requested but not all studies have Nsamples column specified.")


    outfile = args.path_to_res

    with open( outfile, 'w' ) as out:

        ## SHOULD I CORRECT "REF" TO "OTH" AND "ALT" TO "EFF"?
        out.write("\t".join(["#CHR","POS","REF","ALT","SNP", studs[0].name + "_beta", studs[0].name + "_sebeta", studs[0].name + "_pval", studs[0].name + "_af_alt", studs[0].name + "_INFO", studs[0].name + "_Nsamples" ]))

        out.write( ("\t" if len(studs[0].extra_cols) else "") + "\t".join( [studs[0].name + "_" + c for c in studs[0].extra_cols] ) )
        ## align to leftmost STUDY
        for oth in studs[1:len(studs)]:
            out.write( "\t" +  "\t".join( [ oth.name + "_beta", oth.name + "_sebeta", oth.name + "_pval", oth.name + "_af_alt", oth.name + "_INFO", oth.name + "_Nsamples" ] ))
            out.write( ("\t" if len(oth.extra_cols) else "") + "\t".join( [oth.name + "_" + c for c in oth.extra_cols] ) )

            if args.pairwise_with_first:
                for m in methods:
                    out.write("\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_beta\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_sebeta\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_p\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_Ncases\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_Ncontrols\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_Neffective\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_af_alt\t" +
                              studs[0].name + "_" + oth.name + "_" +  m + "_meta_info")

        out.write("\tall_meta_Nstudies\tall_meta_Nsamples")
        for m in methods:
            if args.is_het_test:
                out.write("\tall_"+m+"_meta_beta\t" +
                          "all_"+m+"_meta_sebeta\t" +
                          "all_"+m+"_meta_p\t" +
                          "all_"+m+"_meta_Ncases\t" +
                          "all_"+m+"_meta_Ncontrols\t" +
                          "all_"+m+"_meta_Neffective\t" +
                          "all_"+m+"_meta_af_alt\t" +
                          "all_"+m+"_meta_info\t" +
                          "all_"+m+"_het_p")
            else:
                out.write("\tall_"+m+"_meta_beta\t" +
                          "all_"+m+"_meta_sebeta\t" +
                          "all_"+m+"_meta_p\t" +
                          "all_"+m+"_meta_Ncases\t" +
                          "all_"+m+"_meta_Ncontrols\t" +
                          "all_"+m+"_meta_Neffective\t" +
                          "all_"+m+"_meta_af_alt\t" +
                          "all_"+m+"_meta_info")


        if args.leave_most_sig_out:
            out.write("\tlmso_meta_removed_study_name\tlmso_meta_Nsamples")
            for m in methods:
                if args.is_het_test:
                    out.write("\tlmso_"+m+"_meta_beta\t" +
                              "lmso_"+m+"_meta_sebeta\t" +
                              "lmso_"+m+"_meta_p\t" +
                              "lmso_"+m+"_meta_Ncases\t" +
                              "lmso_"+m+"_meta_Ncontrols\t" +
                              "lmso_"+m+"_meta_Neffective\t" +
                              "lmso_"+m+"_meta_af_alt\t" +
                              "lmso_"+m+"_meta_info\t" +
                              "lmso_"+m+"_het_p")
                else:
                    out.write("\tlmso_"+m+"_meta_beta\t" +
                              "lmso_"+m+"_meta_sebeta\t" +
                              "lmso_"+m+"_meta_p\t" +
                              "lmso_"+m+"_meta_Ncases\t" +
                              "lmso_"+m+"_meta_Ncontrols\t" +
                              "lmso_"+m+"_meta_Neffective\t" +
                              "lmso_"+m+"_meta_af_alt\t" +
                              "lmso_"+m+"_meta_info\t")


        if args.leave_one_out:
            for s in studs:
                out.write("\t" + "leave_" + s.name + "_Nsamples")
                for m in methods:
                    if args.is_het_test:
                        out.write( "\t" +  "\t".join( ["leave_" + s.name + "_" + m + "_meta_beta",
                                                       "leave_" + s.name + "_" + m + "_meta_sebeta",
                                                       "leave_" + s.name + "_" + m + "_meta_p",
                                                       "leave_" + s.name + "_" + m + "_meta_Ncases",
                                                       "leave_" + s.name + "_" + m + "_meta_Ncontrols",
                                                       "leave_" + s.name + "_" + m + "_meta_Neffective",
                                                       "leave_" + s.name + "_" + m + "_meta_af_alt",
                                                       "leave_" + s.name + "_" + m + "_meta_info",
                                                       "leave_" + s.name + "_" + m + "_meta_het_p"] ))
                    else:
                        out.write( "\t" +  "\t".join( ["leave_" + s.name + "_" + m + "_meta_beta",
                                                       "leave_" + s.name + "_" + m + "_meta_sebeta",
                                                       "leave_" + s.name + "_" + m + "_meta_p",
                                                       "leave_" + s.name + "_" + m + "_meta_Ncases",
                                                       "leave_" + s.name + "_" + m + "_meta_Ncontrols",
                                                       "leave_" + s.name + "_" + m + "_meta_Neffective",
                                                       "leave_" + s.name + "_" + m + "_meta_af_alt",
                                                       "leave_" + s.name + "_" + m + "_meta_info"] ))

        out.write("\n")

        next_var = get_next_variant(studs)
        if not args.quiet:
            print("NEXT VARIANTS")
            for v in next_var:
                print(v)
        matching_studies = [(studs[i],v) for i,v in enumerate(next_var) if v is not None]

        while len(matching_studies)>0:

            # Variant info - 5 columns
            d = matching_studies[0][1]
            outdat = [ d.chr, d.pos, d.ref, d.alt]
            v = "{}:{}:{}:{}".format(*outdat)
            outdat.append(v)

            for i,_ in enumerate(studs):
                if next_var[i] is not None:
                    # study stats for variant (if present) - 6 columns + nExtra
                    outdat.extend([format_num(next_var[i].beta), format_num(next_var[i].se), format_num(next_var[i].pval), format_num(next_var[i].af_alt), format_num(next_var[i].info), int(next_var[i].Nsamples) ])
                    outdat.extend([ c for c in next_var[i].extra_cols ])

                    # meta analyse pairwise only with the leftmost study
                    if not args.pairwise_with_first or i==0:
                        continue

                    if next_var[0] is not None:
                        met = do_meta( [(studs[0],next_var[0]), (studs[i],next_var[i])], methods=methods, is_het_test=args.is_het_test)
                        for m in met:
                            if args.is_het_test:
                                # pairwise meta results - 9 columns (per method) if het test
                                outdat.extend([format_num(num) for num in m[0:3]])
                                outdat.extend([int(num) for num in m[3:6]])
                                outdat.extend([format_num(num) for num in m[6:9]])
                            else:
                                # pairwise meta results - 8 columns (per method) if no het test
                                outdat.extend([format_num(num) for num in m[0:3]])
                                outdat.extend([int(num) for num in m[3:6]])
                                outdat.extend([format_num(num) for num in m[6:8]])

                    else:
                        if args.is_het_test:
                            # pairwise not possible: left-most study is missing var: 9 missing columns if het test
                            outdat.extend(["NA"] * len(methods) * 9)
                        else:
                            # pairwise not possible: left-most study is missing var: 8 missing columns if no het test
                            outdat.extend(["NA"] * len(methods) * 8)
                else:
                    # SHOULD THESE BE THE OPPOSITE WAY ROUND? I HAVE CORRECTED FOR NOW, BUT SWITCH 9 & 8 IF NEEDED
                    if args.is_het_test:
                        # pairwise not possible for this study - (6+nExtra)+(9*nMethods) empty cols if pairwise and het test and not first study, otherwise (6+nExtra) empty cols if pairwise and het test and first study
                        outdat.extend(['NA'] * (6 + len(studs[i].extra_cols) + (len(methods)*9 if args.pairwise_with_first and i>0 else 0) ) )
                    else:
                        # pairwise not possible for this study - (6+nExtra)+(9*nMethods) empty cols if pairwise and no het test and not first study, otherwise (6+nExtra) empty cols if pairwise and no het test and first study
                        outdat.extend(['NA'] * (6 + len(studs[i].extra_cols) + (len(methods)*8 if args.pairwise_with_first and i>0 else 0) ) )

            meta_res = []
            if len( matching_studies )>1:
                met = do_meta( matching_studies, methods=methods, is_het_test=args.is_het_test )
                for m in met:
                    if m is not None:
                        if args.is_het_test:
                            # standard meta results - 9 columns (per method) if het test
                            meta_res.extend([format_num(num) for num in m[0:3]])
                            meta_res.extend([int(num) for num in m[3:6]])
                            meta_res.extend([format_num(num) for num in m[6:9]])
                        else:
                            # standard meta results - 8 columns (per method) if no het tes
                            meta_res.extend([format_num(num) for num in m[0:3]])
                            meta_res.extend([int(num) for num in m[3:6]])
                            meta_res.extend([format_num(num) for num in m[6:8]])
                    else:
                        if args.is_het_test:
                            # standard meta returned None so 9 missing columns (per method) if het test
                            meta_res.extend(['NA'] * 9)
                        else:
                            # standard meta returned None so 8 missing columns (per method) if no het test
                            meta_res.extend(['NA'] * 8)


            else:
                if args.is_het_test:
                    # only one study has results - per method: 8 columns (var sumstats) for left-most study and NA column if het test
                    meta_res.extend( [format_num(matching_studies[0][1].beta),
                                      format_num(matching_studies[0][1].se) ,
                                      format_num(matching_studies[0][1].pval),
                                      int(matching_studies[0][0].n_cases) ,
                                      int(matching_studies[0][0].n_controls) ,
                                      int(matching_studies[0][0].effective_size_logistic),
                                      format_num(matching_studies[0][1].af_alt),
                                      format_num(matching_studies[0][1].info),
                                      'NA']  * len(methods) )
                else:
                    # only one study has results - per method: 8 columns (var sumstats) for left-most study if no het test
                    meta_res.extend( [format_num(matching_studies[0][1].beta),
                                      format_num(matching_studies[0][1].se) ,
                                      format_num(matching_studies[0][1].pval),
                                      int(matching_studies[0][0].n_cases) ,
                                      int(matching_studies[0][0].n_controls) ,
                                      int(matching_studies[0][0].effective_size_logistic),
                                      format_num(matching_studies[0][1].af_alt),
                                      format_num(matching_studies[0][1].info)]  * len(methods) )

            # single column - number of studies with the current variant
            outdat.append( str(len(matching_studies)) )
            # single column - number of total samples in meta for current variant
            outdat.append( sum([matching_studies[j][1].Nsamples for j in range(len(matching_studies))]) )
            # append meta-analysis results (or missing columns)
            outdat.extend(meta_res)

            # leave most-significant result out
            if args.leave_most_sig_out:
                # only run if at least two studies left after removing most-sig
                if len( matching_studies )>1:
                    # get number of study for most significant p-value
                    #mspsi = min(range(len(matching_studies)), key = lambda j: matching_studies[j][1].pval)
                    mspsi = numpy.argmin( [matching_studies[j][1].pval for j in range(len(matching_studies))] )
                    # extract matching study results without most-significant study
                    matching_studies_lmso = [(studs[i], var) for i,var in enumerate(next_var) if studs[i].name != matching_studies[mspsi][0].name and var is not None]
                    # add removed study name and sample size for lmso meta column
                    outdat.append( matching_studies[mspsi][0].name )
                    outdat.append( matching_studies_lmso[0][1].Nsamples )

                    # run meta-analysis if 3 or more studies have data for this variant
                    if len( matching_studies )>2:

                        # run meta-analysis with these results
                        met = do_meta( matching_studies_loo, methods=methods, is_het_test=args.is_het_test )
                        # print meta-analysis results for each method
                        for m in met:
                            if m is not None:
                                # if het test flag is set, 9 values outputted from meta
                                if args.is_het_test:
                                    outdat.extend([format_num(num) for num in m[0:3]])
                                    outdat.extend([int(num) for num in m[3:6]])
                                    outdat.extend([format_num(num) for num in m[6:9]])
                                # otherwise, 8 values outputted
                                else:
                                    outdat.extend([format_num(num) for num in m[0:3]])
                                    outdat.extend([int(num) for num in m[3:6]])
                                    outdat.extend([format_num(num) for num in m[6:8]])
                            else:
                                # m returned as None (meta not performed for this method for this var)
                                # add NA*9 if het test flag set, 8*NA otherwise
                                if args.is_het_test:
                                    outdat.extend(['NA'] * 9)
                                else:
                                    outdat.extend(['NA'] * 8)

                    # if only two matching studies, then just print stats of least significant study
                    elif len( matching_studies )==2:

                        # use stats from only study remaining - 9 columns per method if het test flag set
                        if args.is_het_test:
                            outdat.extend( [format_num(matching_studies_lmso[0][1].beta),
                                            format_num(matching_studies_lmso[0][1].se),
                                            format_num(matching_studies_lmso[0][1].pval),
                                            int(matching_studies_lmso[0][0].n_cases),
                                            int(matching_studies_lmso[0][0].n_controls),
                                            int(matching_studies_lmso[0][0].effective_size_logistic) ,
                                            format_num(matching_studies_lmso[0][1].af_alt),
                                            format_num(matching_studies_lmso[0][1].info),
                                            'NA']  * len(methods) )

                        # use stats from only study remaining - 8 columns per method if het test flag not set
                        else:
                            outdat.extend( [format_num(matching_studies_lmso[0][1].beta),
                                            format_num(matching_studies_lmso[0][1].se),
                                            format_num(matching_studies_lmso[0][1].pval),
                                            int(matching_studies_lmso[0][0].n_cases),
                                            int(matching_studies_lmso[0][0].n_controls),
                                            int(matching_studies_lmso[0][0].effective_size_logistic) ,
                                            format_num(matching_studies_lmso[0][1].af_alt),
                                            format_num(matching_studies_lmso[0][1].info)]  * len(methods) )

                # else, only one study so not possible to remove - set lmso columns to NA
                else:
                    # add NA for each removed study name column and lmso Nsamples column
                    outdat.extend(['NA'] * 2)
                    if args.is_het_test:
                        # 9 columns if het test flag set
                        outdat.extend(['NA'] * 9 * len(methods))
                    else:
                        # 8 columns if het test flag not set
                        outdat.extend(['NA'] * 8 * len(methods))


            # leave one out results
            if args.leave_one_out:
                for s,_ in enumerate(studs):
                    matching_studies_loo = [(studs[i], var) for i,var in enumerate(next_var) if s != i and var is not None]
                    #outdat.append( str(len(matching_studies_loo)) )
                    outdat.append( sum([matching_studies_loo[j][1].Nsamples for j in range(len(matching_studies_loo))]) )
                    if len(matching_studies_loo) > 1:
                        met = do_meta( matching_studies_loo, methods=methods, is_het_test=args.is_het_test )
                        for m in met:
                            if m is not None:
                                if args.is_het_test:
                                    outdat.extend([format_num(num) for num in m[0:3]])
                                    outdat.extend([int(num) for num in m[3:6]])
                                    outdat.extend([format_num(num) for num in m[6:9]])
                                else:
                                    outdat.extend([format_num(num) for num in m[0:3]])
                                    outdat.extend([int(num) for num in m[3:6]])
                                    outdat.extend([format_num(num) for num in m[6:8]])
                            else:
                                if args.is_het_test:
                                    outdat.extend(['NA'] * 9)
                                else:
                                    outdat.extend(['NA'] * 8)

                    elif len(matching_studies_loo) == 1:
                        if args.is_het_test:
                            outdat.extend( [format_num(matching_studies_loo[0][1].beta),
                                            format_num(matching_studies_loo[0][1].se) ,
                                            format_num(matching_studies_loo[0][1].pval),
                                            int(matching_studies_loo[0][0].n_cases) ,
                                            int(matching_studies_loo[0][0].n_controls) ,
                                            int(matching_studies_loo[0][0].effective_size_logistic) ,
                                            format_num(matching_studies_loo[0][1].af_alt),
                                            format_num(matching_studies_loo[0][1].info),
                                            'NA']  * len(methods) )
                        else:
                            outdat.extend( [format_num(matching_studies_loo[0][1].beta),
                                            format_num(matching_studies_loo[0][1].se) ,
                                            format_num(matching_studies_loo[0][1].pval),
                                            int(matching_studies_loo[0][0].n_cases) ,
                                            int(matching_studies_loo[0][0].n_controls) ,
                                            int(matching_studies_loo[0][0].effective_size_logistic),
                                            format_num(matching_studies_loo[0][1].af_alt),
                                            format_num(matching_studies_loo[0][1].info)]  * len(methods) )
                    else:
                        if args.is_het_test:
                            outdat.extend(['NA'] * 9 * len(methods))
                        else:
                            outdat.extend(['NA'] * 8 * len(methods))

            out.write( "\t".join([ str(o) for o in outdat]) + "\n" )

            next_var = get_next_variant(studs)
            if not args.quiet:
                print("NEXT VARIANTS")
                for v in next_var:
                    print(v)
            matching_studies = [(studs[i],v) for i,v in enumerate(next_var) if v is not None]

    subprocess.run(["bgzip","--force",args.path_to_res])
    subprocess.run(["tabix","-s 1","-b 2","-e 2",args.path_to_res + ".gz"])


if __name__ == '__main__':
    run()
