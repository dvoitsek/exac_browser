import re
from utils import *
import log
import bson

SEARCH_LIMIT = 10000


def get_gene(db, gene_id):
    return db.genes.find_one({'gene_id': gene_id}, projection={'_id': False})


def get_gene_by_name(db, gene_name):
    # try gene_name field first
    gene = db.genes.find_one({'gene_name': gene_name}, projection={'_id': False})
    if gene:
        return gene
    # if not, try gene['other_names']
    return db.genes.find_one({'other_names': gene_name}, projection={'_id': False})


def get_transcript(db, transcript_id):
    transcript = db.transcripts.find_one({'transcript_id': transcript_id}, projection={'_id': False})
    if not transcript:
        return None
    transcript['exons'] = get_exons_in_transcript(db, transcript_id)
    return transcript


def get_raw_variant(db, xpos, ref, alt, get_id=False):
    # log.print_info("Search with xpos=%d ref=%s alt=%s" % (xpos, ref, alt))
    variant = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection={'_id': get_id})
    if variant is not None:
        expand_variant_annotations(variant, db)
    else:
        log.print_warning("Found no variant for xpos=%d ref=%s alt=%s" % (xpos, ref, alt))
    return variant


def get_variant(db, xpos, ref, alt):
    variant = get_raw_variant(db, xpos, ref, alt, False)
    if variant is None or 'rsid' not in variant:
        return variant
    if variant['rsid'] == '.' or variant['rsid'] is None:
        rsid = db.dbsnp.find_one({'xpos': xpos})
        if rsid:
            variant['rsid'] = 'rs%s' % rsid['rsid']
    return variant


def get_variants_by_rsid(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    variants = list(db.variants.find({'rsid': rsid}, projection={'_id': False}))
    variants = [expand_variant_annotations(variant, db) for variant in variants]
    add_consequence_to_variants(variants)
    variants = [pick_variant_source(variant, -1) for variant in variants]
    return variants


def get_variants_from_dbsnp(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        rsid = int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    position = db.dbsnp.find_one({'rsid': rsid})
    if position:
        variants = list(db.variants.find({'xpos': {'$lte': position['xpos'], '$gte': position['xpos']}}, projection={'_id': False}))
        if variants:
            variants = [expand_variant_annotations(variant, db) for variant in variants]
            variants = [pick_variant_source(variant, -1) for variant in variants]
            add_consequence_to_variants(variants)
            return variants
    return []


def get_coverage_for_bases(db, xstart, xstop=None):
    """
    Get the coverage for the list of bases given by xstart->xstop, inclusive
    Returns list of coverage dicts
    xstop can be None if just one base, but you'll still get back a list
    """
    if xstop is None:
        xstop = xstart
    coverages = {
        doc['xpos']: doc for doc in db.base_coverage.find(
            {'xpos': {'$gte': xstart, '$lte': xstop}},
            projection={'_id': False}
        )
    }
    ret = []
    for i in range(xstart, xstop+1):
        if i in coverages:
            ret.append(coverages[i])
        else:
            ret.append({'xpos': i, 'pos': xpos_to_pos(i)})
    for item in ret:
        item['has_coverage'] = 'mean' in item
        del item['xpos']
    return ret


def get_coverage_for_transcript(db, xstart, xstop=None):
    """

    :param db:
    :param genomic_coord_to_exon:
    :param xstart:
    :param xstop:
    :return:
    """
    coverage_array = get_coverage_for_bases(db, xstart, xstop)
    # only return coverages that have coverage (if that makes any sense?)
    # return coverage_array
    covered = [c for c in coverage_array if c['has_coverage']]
    for c in covered:
        del c['has_coverage']
    return covered


def get_constraint_for_transcript(db, transcript):
    return db.constraint.find_one({'transcript': transcript}, projection={'_id': False})


def get_exons_cnvs(db, transcript_name):
   return list(db.cnvs.find({'transcript': transcript_name}, projection={'_id': False}))

def get_cnvs(db, gene_name):
   return list(db.cnvgenes.find({'gene': gene_name}, projection={'_id': False}))


def get_awesomebar_suggestions(g, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = [r for r in g.autocomplete_strings if regex.match(r)][:20]
    return results


# 1:1-1000
R1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
R2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
R3 = re.compile(r'^(\d+|X|Y|M|MT)$')
# R4 = re.compile(r'^(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)-([ATCG]+)-([ATCG]+)$')
R4 = re.compile(r'^\s*(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)[-:\s]*([ATCG]+)\s*[-:/]\s*([ATCG]+)\s*$')


def get_awesomebar_result(db, query):
    """
    Similar to the above, but this is after a user types enter
    We need to figure out what they meant - could be gene, variant, region

    Return tuple of (datatype, identifier)
    Where datatype is one of 'gene', 'variant', or 'region'
    And identifier is one of:
    - ensembl ID for gene
    - variant ID string for variant (eg. 1-1000-A-T)
    - region ID string for region (eg. 1-1000-2000)

    Follow these steps:
    - if query is an ensembl ID, return it
    - if a gene symbol, return that gene's ensembl ID
    - if an RSID, return that variant's string


    Finally, note that we don't return the whole object here - only it's identifier.
    This could be important for performance later

    """
    query = query.strip()
    print 'Query: %s' % query

    # Variant
    variant = get_variants_by_rsid(db, query.lower())
    if variant:
        if len(variant) == 1:
            return 'variant', variant[0]['variant_id']
        else:
            return 'dbsnp_variant_set', variant[0]['rsid']
    variant = get_variants_from_dbsnp(db, query.lower())
    if variant:
        return 'variant', variant[0]['variant_id']
    # variant = get_variant(db, )
    # TODO - https://github.com/brettpthomas/exac_browser/issues/14

    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']

    # From here out, all should be uppercase (gene, tx, region, variant_id)
    query = query.upper()
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']

    # Ensembl formatted queries
    if query.startswith('ENS'):
        # Gene
        gene = get_gene(db, query)
        if gene:
            return 'gene', gene['gene_id']

        # Transcript
        transcript = get_transcript(db, query)
        if transcript:
            return 'transcript', transcript['transcript_id']

    # From here on out, only region queries
    if query.startswith('CHR'):
        query = query.lstrip('CHR')
    # Region
    m = R1.match(query)
    if m:
        if int(m.group(3)) < int(m.group(2)):
            return 'region', 'invalid'
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))
    m = R2.match(query)
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(2))
    m = R3.match(query)
    if m:
        return 'region', '{}'.format(m.group(1))
    m = R4.match(query)
    if m:
        return 'variant', '{}-{}-{}-{}'.format(m.group(1), m.group(2), m.group(3), m.group(4))

    return 'not_found', query


def get_genes_in_region(db, chrom, start, stop):
    """
    Genes that overlap a region
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    genes = db.genes.find({
        'xstart': {'$lte': xstop},
        'xstop': {'$gte': xstart},
    }, projection={'_id': False})
    return list(genes)


def get_variants_in_region(db, chrom, start, stop):
    """
    Variants that overlap a region
    Unclear if this will include CNVs
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    variants = list(db.variants.find({
        'xpos': {'$lte': xstop, '$gte': xstart}
    }, projection={'_id': False}, limit=SEARCH_LIMIT))
    variants = [expand_variant_annotations(variant, db) for variant in variants]
    variants = [pick_variant_source(variant, -1) for variant in variants]
    add_consequence_to_variants(variants)
    for variant in variants:
        # remove_extraneous_information(variant)
        strip_internal_fields(variant)
    return list(variants)


def get_metrics(db, variant):
    if 'allele_count' not in variant or variant['allele_num'] == 0:
        return None
    metrics = {}
    for metric in METRICS:
        metrics[metric] = db.metrics.find_one({'metric': metric}, projection={'_id': False})

    metric = None
    if variant['allele_count'] == 1:
        metric = 'singleton'
    elif variant['allele_count'] == 2:
        metric = 'doubleton'
    else:
        for af in AF_BUCKETS:
            if float(variant['allele_count'])/variant['allele_num'] < af:
                metric = af
                break
    if metric is not None:
        metrics['Site Quality'] = db.metrics.find_one({'metric': 'binned_%s' % metric}, projection={'_id': False})
    return metrics


def remove_extraneous_information(variant):
    del variant['genotype_depths']
    del variant['genotype_qualities']
    del variant['transcripts']
    del variant['genes']
    del variant['orig_alt_alleles']
    del variant['xpos']
    del variant['xstart']
    del variant['xstop']
    del variant['site_quality']
    del variant['vep_annotations']


def get_variants_in_gene(db, gene_id):
    """
    """
    variants = []
    annotation_ids = []
    for annotation in db.annotations.find({'$or' : [{'Gene_Name': gene_id}, {'Gene_ID': gene_id}]}):
        annotation_ids.append(str(annotation['_id']))
    # log.print_info("Found %d annotations" % len(annotation_ids))
    for variant in db.variants.find({'ANN._id': {'$in' : annotation_ids}}, projection={'_id': False}):
        # variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Gene'] == gene_id]
        pick_variant_source(variant, -1)
        expand_variant_annotations(variant, db)
        add_consequence_to_variant(variant)
        strip_internal_fields(variant)
        # remove_extraneous_information(variant)
        variants.append(variant)
    # log.print_info("Found %d variants matching gene id or name" % len(variants))
    return variants


def get_transcripts_in_gene(db, gene_id):
    """
    """
    return list(db.transcripts.find({'gene_id': gene_id}, projection={'_id': False}))


def get_variants_in_transcript(db, transcript_id):
    """
    """
    variants = []
    annotation_ids = []
    # Get transcript's gene ids
    gene_ids = list(db.transcripts.find({'transcript_id': transcript_id}, projection={'_id': False, 'gene_id': True}))
    gene_ids = [x['gene_id'] for x in gene_ids]
    log.print_info("Found genes in transcript %s" % gene_ids)

    # Find annotations mentioning gene ids
    for annotation in db.annotations.find({'$or' : [{'Gene_Name': {'$in': gene_ids}}, {'Gene_ID': {'$in': gene_ids}}]}):
        annotation_ids.append(str(annotation['_id']))

    # Grab corresponding variants
    for variant in db.variants.find({'ANN._id': {'$in' : annotation_ids}}, projection={'_id': False}):
        # variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Feature'] == transcript_id]
        pick_variant_source(variant, -1)
        expand_variant_annotations(variant, db)
        add_consequence_to_variant(variant)
        strip_internal_fields(variant)
        # remove_extraneous_information(variant)
        variants.append(variant)
    return variants


def get_exons_in_transcript(db, transcript_id):
    # return sorted(
    #     [x for x in
    #      db.exons.find({'transcript_id': transcript_id}, fields={'_id': False})
    #      if x['feature_type'] != 'exon'],
    #     key=lambda k: k['start'])
    return sorted(list(db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, projection={'_id': False})), key=lambda k: k['start'])

def expand_variant_annotations(variant, db):
    '''
    Input variant MUST be unsourced
    '''
    annotations = []
    for ann in variant['ANN']:
        # log.print_info("Expanding annotation %s" % ann)
        src_ann = db.annotations.find_one(bson.objectid.ObjectId(oid=ann['_id']), projection={'_id': False})
        src_ann['_source'] = ann['_source']
        if src_ann is None:
            log.print_error("Annotation %s referenced in variant %s not found in annotations collection." % (ann['_id'], variant['variant_id']))
        else:
            annotations.append(src_ann)
    variant['ANN'] = annotations
    return variant
