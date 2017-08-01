#!/usr/bin/env python
'''Help for running this stuff.
How does it work?
'''

# from flask.ext.script import Manager
from flask_script import Manager
from exac import app
import exac


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
@manager.option('-d', '--directory', dest='directory', default=None)
def load_base_coverage(directory):
    return if directory == None
    exac.load_base_coverage(directory)


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def load_variants_file(file):
    return if file == None
    exac.load_variants_file(file)


# @manager.command
# def reload_variants():
#     exac.load_variants_file()
#     exac.load_mnps()
#     exac.precalculate_metrics()


@manager.command
@manager.option('-f', '--file', dest='file', default=None)
def load_gene_models(file):
    return if file == None
    exac.load_gene_models(file)


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
@manager.option('-f', '--file', dest='file', default=None)
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

