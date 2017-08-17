#!/usr/bin/env python
'''Help for running this stuff.
How does it work?
'''

# from flask.ext.script import Manager
from flask_script import Manager
from exac import app
import exac


ACCEPTABLE_REF_ALTS = ['GRCh37-']
manager = Manager(app)


@manager.command
def hello():
    '''
    Hello world
    '''
    print "hello"


# @manager.command
# def load_db():
#    exac.load_db()


@manager.command
@manager.option('-d', '--directory', dest='directory', default=None, help="Directory containing coverage files. All files in the directory will be loaded. Please make sure it contains only coverage files.")
@manager.option('-f', '--files', dest='files', default=[], help="List of coverage files to load.")
def load_base_coverage(directory, *files):
    '''
    Loads coverage files specified by the input.
    If directory is provided, the script will attempt to load all files contained in the directory.
    If both directory and files are provided the script will attempt to load both the content of the directory and of the files.
    All input must be gzipped.
    '''
    raise IOError("No input provided.") if directory == None && len(files) == 0
    exac.load_base_coverage(directory, files)


@manager.command
@manager.option('-r', '--reference', dest='reference', required=True)
@manager.option('-a', '--reference-alternative', dest='alt', required=False, default='')
@manager.option('-f', '--files', dest='files', default=[], required=True)
def load_variants_file(reference, alt, *files):
    '''
    Loads VCF files specified by the input. All input must be gzipped and tabix-indexed.
    Currently accepted input references and alternatives:
        - GRCh37
    If your reference does not appear in this list, contact voitsekh@cng.fr
    '''
    ref_alt_string = "%s-%s"(reference, alt)
    raise IOError("No input provided.") if len(files) == 0
    raise IOError("Not acceptable reference-alternative combination %s"(ref_alt_string)) if ref_alt_string not in ACCEPTABLE_REF_ALTS
    exac.load_variants_file(reference, alt, files)

@manager.command
@manager.option('-f', '--files', dest='files', default=[], required=True)
def drop_sources(*files):
    '''
    Drop all data associated to this file from the database.
    '''
    raise IOError('No input provided') if len(files) == 0
    exac.drop_sources(files)

# @manager.command
# def reload_variants():
#     exac.load_variants_file()
#     exac.load_mnps()
#     exac.precalculate_metrics()


@manager.command
@manager.option('-c', '--canonical_transcripts', dest='canonical_transcripts', default=None, required=False)
@manager.option('-o', '--omim', dest='omim', default=None, required=False)
@manager.option('-d', '--dbnsfp', dest='dbnsfp', default=None, required=False)
@manager.option('-g', '--gencode', dest='gencode', default=None, required=False)
def load_gene_models(canonical_transcripts, omim, dbnsfp, gencode):
    """
    Loads gene models from 4 files.
    """
    exac.load_gene_models(canonical_transcripts, omim, dbnsfp, gencode)


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def load_cnv_models(file):
    return if file == None
    exac.load_cnv_models(file)


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def load_cnv_genes(file):
    return if file == None
    exac.load_cnv_genes(file)


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def drop_cnv_genes(file):
    return if file == None
    exac.drop_cnv_genes(file)


@manager.command
@manager.option('-f', '--file', dest='file', default=None, required=True)
def load_dbsnp_file(file):
    return if file == None
    exac.load_dbsnp_file(file)


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def load_constraint_information(file):
    return if file == None
    exac.load_constraint_information(file)


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def load_mnps(file):
    return if file == None
    exac.load_mnps(file)


@manager.command
def create_cache():
    exac.create_cache()


@manager.command
def precalculate_metrics():
    exac.precalculate_metrics()

if __name__ == "__main__":
    manager.run()
