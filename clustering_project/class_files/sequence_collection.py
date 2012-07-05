#!/usr/bin/python
# -*- coding: utf-8 -*-

import glob
import os
import multiprocessing
import copy
import re
import matplotlib.pyplot as plt
from sequence_record import TCSeqRec
from clustering import Clustering
import copy_reg
import types
from random import shuffle
import shutil


def _pickle_method(method):
    """
    Adjust pickling via copy_reg module to make multiprocessing.Pool work
    with class methods (otherwise unpickleable)
    """

    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):

        # deal with mangled names

        cls_name = cls.__name__.lstrip('_')
        func_name = '_%s%s' % (cls_name, func_name)
    return (_unpickle_method, (func_name, obj, cls))


def _unpickle_method(func_name, obj, cls):
    """
    Adjust pickling via copy_reg module to make multiprocessing.Pool work
    with class methods (otherwise unpickleable)
    """

    if obj and func_name in obj.__dict__:
        (cls, obj) = (obj, None)  # if func_name is classmethod
    for cls in cls.__mro__:
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


class SequenceCollection(object):

    """
    Orchestrating class that should:
    a) work as a central repository for the information generated by the 
       subordinate classes, and
    b) be the only class directly interacted with by the user 

    TO DO:
    implement consistent naming of methods (where appropriate)
    Prefixes:
    get_[something]  - returns the object implied by something
    put_[something]  - puts something in the class data structure
    show_[something] - prints something to screen
    plot_[something] - displays a plot of something
    _[something]     - private method
    """

    def get_files(self, input_dir, file_format='fasta'):
        """
        Get list of alignment files from an input directory
        *.fa, *.fas and *.phy files only
        Stores in self.files
        """

        if file_format == 'fasta':
            files = glob.glob('{0}/*.fa'.format(input_dir))
            if len(files) == 0:
                files = glob.glob('{0}/*.fas'.format(input_dir))
        elif file_format == 'phylip':
            files = glob.glob('{0}/*.phy'.format(input_dir))
        else:
            print 'Unrecognised file format %s' % file_format
            files = None
        if not files:
            print 'No sequence files found in {0}'.format(input_dir)
            return 0
        return sorted(files)

    def put_records(
        self,
        files,
        file_format='fasta',
        datatype='protein',
        ):
        """ 
        Reads sequence files from the list generated by
        get_files and stores in self.records
        """

        name = lambda i: i[i.rindex('/') + 1:i.rindex('.')]
        self.records = [TCSeqRec(f, file_format=file_format,
                        name=name(f), datatype=datatype) for f in files]

    def get_records(self):
        """
        Returns list of stored sequence records
        """

        return self.records

    def sanitise_records(self):
        for rec in self.records:
            rec.sanitise()

    def put_dv_matrices(
        self,
        tmpdir='/tmp',
        helper='./class_files/DV_wrapper.drw',
        overwrite=True,
        ):

        for rec in self.records:
            rec.dv = [rec.get_dv_matrix(tmpdir=tmpdir, helper=helper,
                      overwrite=overwrite)]

    def _unpack_dv(self, packed_args):
        return packed_args[0].get_dv_matrix(*packed_args[1:])

    def _dv_parallel_call(
        self,
        tmpdir='/tmp',
        helper='./class_files/DV_wrapper.drw',
        overwrite=True,
        ):

        nprocesses = min(len(self.records), multiprocessing.cpu_count()
                         - 1)
        print 'Initialising a pool of {0} processes running {1} jobs...'.format(nprocesses,
                len(self.records))
        pool = multiprocessing.Pool(nprocesses)
        results = []
        args = []
        names = []
        for rec in self.records:
            new_dir = tmpdir + '/' + rec.name
            if not os.path.isdir(new_dir):
                os.mkdir(new_dir)
            args.append((rec, tmpdir + '/' + rec.name, helper,
                        overwrite))
            names.append(rec.name)
        r = pool.map_async(self._unpack_dv, args,
                           callback=results.append)
        r.wait()
        for (w, x, y, z) in args:
            if os.path.isdir(x):
                os.rmdir(x)
        results = results[0]
        print 'Results obtained, closing pool...'
        pool.close()
        pool.join()
        print 'Pool closed'
        return dict(zip(names, results))

    def put_dv_matrices_parallel(
        self,
        tmpdir='/tmp',
        helper='./class_files/DV_wrapper.drw',
        overwrite=True,
        ):

        dv_matrices_dict = self._dv_parallel_call(tmpdir, helper,
                overwrite=overwrite)
        for rec in self.records:
            rec.dv = [dv_matrices_dict[rec.name]]

    def get_dv_matrices(self):
        dvs = {}
        for rec in self.get_records():
            dvs[rec.name] = rec.dv
        return dvs

    def _unpack_phyml(self, packed_args):
        return packed_args[0].get_phyml_tree(*packed_args[1:])

    def _phyml_parallel_call(
        self,
        model=None,
        datatype=None,
        rec_list=None,
        ncat=4,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if not rec_list:
            rec_list = self.records
        nprocesses = min(len(rec_list), multiprocessing.cpu_count() - 1)
        print 'Initialising a pool of {0} processes running {1} jobs...'.format(nprocesses,
                len(rec_list))
        pool = multiprocessing.Pool(nprocesses)
        results = []
        args = []
        names = []
        for rec in rec_list:
            args.append((
                rec,
                model,
                datatype,
                ncat,
                tmpdir,
                overwrite,
                ))
            names.append(rec.name)
        r = pool.map_async(self._unpack_phyml, args,
                           callback=results.append)
        r.wait()
        print 'Results obtained, closing pool...'
        pool.close()
        pool.join()
        print 'Pool closed'
        return dict(zip(names, results[0]))

    def _unpack_raxml(self, packed_args):
        return packed_args[0].get_raxml_tree(*packed_args[1:])

    def _raxml_parallel_call(
        self,
        rec_list=None,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if not rec_list:
            rec_list = self.records
        nprocesses = multiprocessing.cpu_count() - 1
        pool = multiprocessing.Pool(nprocesses)
        results = []
        args = []
        names = []
        for rec in rec_list:
            args.append((rec, tmpdir, overwrite))
            names.append(rec.name)
        r = pool.map_async(self._unpack_raxml, args,
                           callback=results.append)
        r.wait()
        pool.close()
        pool.join()
        return dict(zip(names, results[0]))

    def _unpack_TC(self, packed_args):
        return packed_args[0].get_TC_tree(*packed_args[1:])

    def _TC_parallel_call(
        self,
        rec_list=None,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if not rec_list:
            rec_list = self.records
        nprocesses = multiprocessing.cpu_count() - 1
        pool = multiprocessing.Pool(nprocesses)
        results = []
        args = []
        names = []
        for rec in rec_list:
            args.append((rec, tmpdir, overwrite))
            names.append(rec.name)
        r = pool.map_async(self._unpack_TC, args,
                           callback=results.append)
        r.wait()
        pool.close()
        pool.join()
        return dict(zip(names, results[0]))

    def put_trees(
        self,
        rec_list=None,
        program='treecollection',
        model=None,
        datatype=None,
        ncat=4,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if not program in ['treecollection', 'raxml', 'phyml']:
            print 'unrecognised program {0}'.format(program)
            return
        if not rec_list:
            rec_list = self.records
        for rec in rec_list:
            if program == 'treecollection':
                rec.get_TC_tree(tmpdir=tmpdir, overwrite=overwrite)
            elif program == 'raxml':
                rec.get_raxml_tree(tmpdir=tmpdir, overwrite=overwrite)
            elif program == 'phyml':
                rec.get_phyml_tree(model=model, datatype=datatype,
                                   tmpdir=tmpdir, ncat=ncat)

    def put_trees_parallel(
        self,
        rec_list=None,
        program='treecollection',
        model=None,
        datatype=None,
        ncat=4,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if not program in ['treecollection', 'raxml', 'phyml']:
            print 'unrecognised program {0}'.format(program)
            return
        if not rec_list:
            rec_list = self.records
        if program == 'treecollection':
            trees_dict = self._TC_parallel_call(rec_list=rec_list,
                    tmpdir=tmpdir, overwrite=overwrite)
        elif program == 'raxml':
            trees_dict = self._raxml_parallel_call(rec_list=rec_list,
                    tmpdir=tmpdir, overwrite=overwrite)
        elif program == 'phyml':
            trees_dict = self._phyml_parallel_call(
                rec_list=rec_list,
                model=model,
                datatype=datatype,
                tmpdir=tmpdir,
                ncat=ncat,
                overwrite=overwrite,
                )
        for rec in self.records:
            rec.tree = trees_dict[rec.name]

    def get_trees(self):
        trees = {}
        for rec in self.records:
            trees[rec.name] = rec.tree
        return trees

    def put_distance_matrices(self, metrics, **kwargs):
        """
        Pass this function a list of metrics
        valid kwargs - invert (bool), normalise (bool)
        """

        if not isinstance(metrics, list):
            metrics = [metrics]
        trees = [rec.tree for rec in self.records]
        for metric in metrics:
            self.clustering.put_distance_matrix(trees, metric, **kwargs)

    def get_distance_matrices(self):
        return self.clustering.distance_matrices

    def put_partitions(
        self,
        metrics,
        linkages,
        nclasses,
        criterion='distance',
        prune=True
        ):
        """
        metrics, linkages and nclasses are given as lists, or coerced into 
        lists
        """

        if not isinstance(metrics, list):
            metrics = [metrics]
        if not isinstance(linkages, list):
            linkages = [linkages]
        if not isinstance(nclasses, list):
            nclasses = [nclasses]
        else:
            nclasses = sorted(nclasses, reverse=True)
        names = [rec.name for rec in self.records]
        for metric in metrics:
            if not metric in self.get_distance_matrices():
                trees = [rec.tree for rec in self.records]  # preserve ordering
                self.clustering.put_distance_matrix(trees, metric)
            for linkage in linkages:
                for n in nclasses:
                    self.clustering.put_partition(metric, linkage, n,
                            names, criterion=criterion, prune=prune)
                    key = (metric, linkage, n)

    def get_partitions(self):
        return self.clustering.partitions

    def put_clusters(self):
        for key in self.clustering.partitions:
            if key not in self.clustering.clusters:
                self.clustering.concatenate_records(key, self.records)

    def get_clusters(self):
        return self.clustering.clusters

    def get_cluster_records(self):
        """
        Returns all concatenated records from cluster analysis
        """

        rec_list = []
        clusters_dict = self.get_clusters()  # This is the outer dictionary

        for compound_key in clusters_dict:
            result = clusters_dict[compound_key]  # Result object holds info
            for i in range(result.length):  # Iterate over clusters in Result object
                rec_list.append(result.retrieve_concat(i))  # Retrieve the concat

        return rec_list

    def put_cluster_trees(
        self,
        program='treecollection',
        model=None,
        datatype=None,
        ncat=4,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if program not in ['treecollection', 'raxml', 'phyml']:
            print 'unrecognised program {0}'.format(program)
            return
        rec_list = self.get_cluster_records()
        self.put_trees(
            rec_list=rec_list,
            program=program,
            model=model,
            ncat=ncat,
            datatype=datatype,
            tmpdir=tmpdir,
            overwrite=overwrite,
            )
        self.update_results()

    def update_results(self):
        for result in self.get_clusters().values():
            result.update_score()

    def put_cluster_trees_parallel(
        self,
        program='treecollection',
        model=None,
        datatype=None,
        ncat=4,
        tmpdir='/tmp',
        overwrite=True,
        ):

        if program not in ['treecollection', 'raxml', 'phyml']:
            print 'unrecognised program {0}'.format(program)
            return
        rec_list = self.get_cluster_records()
        if program == 'treecollection':
            cluster_trees_dict = \
                self._TC_parallel_call(rec_list=rec_list,
                    tmpdir=tmpdir, overwrite=overwrite)
        elif program == 'raxml':
            cluster_trees_dict = \
                self._raxml_parallel_call(rec_list=rec_list,
                    tmpdir=tmpdir, overwrite=overwrite)
        elif program == 'phyml':
            cluster_trees_dict = self._phyml_parallel_call(
                rec_list=rec_list,
                model=model,
                datatype=datatype,
                ncat=ncat,
                tmpdir=tmpdir,
                overwrite=overwrite,
                )
        for rec in rec_list:
            rec.tree = cluster_trees_dict[rec.name]
        self.update_results()

    def get_cluster_trees(self):
        rec_list = sorted(self.get_cluster_records(), key=lambda rec: \
                          rec.name)
        trees = [rec.tree for rec in rec_list]
        return trees

    def get_randomised_alignments(self):

        def pivot(list):
            length = len(list[0])
            new_list = []
            for i in range(length):
                new_list.append(''.join([x[i] for x in list]))
            return new_list

        lengths = [rec.seqlength for rec in self.records]
        datatype = self.records[0].datatype
        concat = copy.deepcopy(self.records[0])
        for rec in self.records[1:]:
            concat += rec
        columns = pivot(concat.sequences)
        shuffle(columns)
        newcols = []
        for l in lengths:
            newcols.append(columns[:l])
            columns = columns[l:]
        newrecs = []
        for col in newcols:
            newseqs = pivot(col)
            newrec = TCSeqRec(headers=concat.headers,
                              sequences=newseqs, datatype=datatype)
            newrecs.append(newrec)
        for i in range(self.length):
            newrecs[i].name = self.records[i].name
        return newrecs

    def make_randomised_copy(self, program='treecollection',
                             tmpdir='/tmp'):
        """
        Lesson learned: don't forget to call pool.close and pool.join
        when multithreading
        """

        other = copy.deepcopy(self)
        other.records = self.get_randomised_alignments()
        other.put_dv_matrices_parallel()
        other.put_trees_parallel(program=program, tmpdir=tmpdir)
        if self.get_clusters():
            other.clustering.clusters = {}
            other.put_clusters()
        if self.get_cluster_trees()[0].newick:
            cluster_program = \
                self.get_cluster_trees()[0].program.lower()
            other.put_cluster_trees_parallel(program=cluster_program,
                    tmpdir=tmpdir)
        return other

    def show_memberships(self):

        partitions = self.get_partitions()
        for compound_key in partitions:
            print ' '.join(str(x) for x in compound_key)
            partition = partitions[compound_key]
            print partition
            print self.clustering.get_memberships(partition)

    def plot_dendrogram(
        self,
        metric,
        link,
        nclasses,
        show=True,
        ):

        plot_object = self.clustering.plot_dendrogram((metric, link,
                nclasses))
        if show:
            plot_object.show()
        return plot_object

    def scree_plot(
        self,
        metric,
        link,
        *args
        ):

        cl = self.get_clusters()
        (x, y) = ([], [])
        for k in sorted(cl):
            if metric in k and link in k:
                y.append(cl[k].score)
                x.append(k[-1])
        plt.plot(x, y, *args)

    def find_mergeable_groups(self, compound_key):

        result_object = self.get_clusters()[compound_key]
        cluster_trees = [rec.tree for rec in result_object.concats]
        cluster_names = [tree.name for tree in cluster_trees]
        matrix = self.clustering.get_distance_matrix(cluster_trees,
                'sym')
        groups = result_object.find_mergeable_groups(matrix)
        return groups

    def merge_groups(self, compound_key):

        group_dict = self.find_mergeable_groups(compound_key)[0]
        old_memberships = self.get_clusters()[compound_key].members
        new_memberships = []
        skips = []
        for i in range(len(old_memberships)):
            if i in skips:
                continue
            elif i in group_dict:
                new = []
                new += old_memberships[i]
                for val in group_dict[i]:
                    new += old_memberships[val]
                    skips.append(val)
                new_memberships.append(new)
                continue
            new_memberships.append(old_memberships[i])

        new_partition = [0 for rec in self.records]
        i = 1
        for group in new_memberships:
            for member in group:
                index = self.records.index(member)
                new_partition[index] = i
            i += 1

        new_key = compound_key + ('merge', )
        self.clustering.partitions[new_key] = new_partition
        self.put_clusters()

    def simulate_from_result(
        self,
        compound_key,
        helper='./class_files/DV_wrapper.drw',
        tmpdir='/tmp',
        ):

        shorten = lambda x: '_'.join([str(b)[:5] for b in x])
        result_object = self.get_clusters()[compound_key]
        lengths = [[rec.seqlength for rec in m] for m in
                   result_object.members]
        total_lengths = [sum(x) for x in lengths]
        msa_dir = '{0}/msa'.format(tmpdir)
        if not os.path.isdir(msa_dir):
            os.mkdir(msa_dir)
        k = 1
        for i in range(result_object.length):
            tree = result_object.concats[i].tree
            tree = tree.pam2sps('sps2pam')
            treefile = \
                tree.write_to_file('{0}/{1}_tmptree{2}.nwk'.format(tmpdir,
                                   shorten(compound_key), i))
            outfile = '{0}/{1}_class{2}_params.drw'.format(tmpdir,
                    shorten(compound_key), i)
            length_list = lengths[i]
            total_length = (total_lengths[i] + total_lengths[i] % 3) / 3
            result_object.write_ALF_parameters(
                'alfsim_{0}'.format(i),
                tmpdir,
                'alftmp_{0}'.format(i),
                1,
                total_length,
                treefile,
                outfile,
                )
            os.system('alfsim {0}'.format(outfile))
            record = \
                TCSeqRec(glob.glob('{0}/alftmp_{1}/alfsim_{1}/MSA/*dna.fa'.format(tmpdir,
                         i))[0])

            alf_newick = \
                open('{0}/alftmp_{1}/alfsim_{1}/RealTree.nwk'.format(tmpdir,
                     i)).read()
            replacement_dict = dict(zip(re.findall(r'(\w+)(?=:)',
                                    alf_newick),
                                    re.findall(r'(\w+)(?=:)',
                                    tree.newick)))
            print alf_newick
            print tree.newick
            print replacement_dict
            record.sort_by_name()
            headers = [replacement_dict[x[:x.rindex('/')]] for x in
                       record.headers]
            print headers
            sequences = record.sequences
            print record
            for j in range(len(length_list)):
                start = sum(length_list[:j])
                end = sum(length_list[:j + 1])
                new_sequences = [seq[start:end] for seq in sequences]
                newmsa = TCSeqRec(headers=headers,
                                  sequences=new_sequences,
                                  name='gene{0:0>3}'.format(k))
                k += 1
                newmsa.write_fasta(outfile='{0}/{1}.fas'.format(msa_dir,
                                   newmsa.name))
            shutil.rmtree('{0}/alftmp_{1}'.format(tmpdir, i))
            os.remove(treefile)
            os.remove(outfile)

        new_seqcol_object = SequenceCollection(msa_dir, datatype='dna',
                helper=helper, tmpdir=tmpdir)
        shutil.rmtree('{0}/msa'.format(tmpdir))
        return new_seqcol_object

    def __init__(
        self,
        input_dir=None,
        records=None,
        file_format='fasta',
        datatype='protein',
        helper='./class_files/DV_wrapper.drw',
        tmpdir='/tmp',
        parallel_load=True,
        overwrite=True,
        ):

        self.dir = input_dir
        self.files = None
        self.records = records
        self.clustering = Clustering()
        self.length = 0

        if input_dir:

            files = self.get_files(input_dir, file_format)
            if files == 0:
                print 'There was a problem reading files from {0}'.format(input_dir)
                return
            if not os.path.isfile(helper):
                print 'There was a problem finding the darwin helper at {0}'.format(helper)
                return
            self.put_records(files, file_format, datatype)
            self.sanitise_records()
            if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir)
            if parallel_load:
                self.put_dv_matrices_parallel(helper=helper,
                        tmpdir=tmpdir, overwrite=overwrite)
            else:
                self.put_dv_matrices(helper=helper, tmpdir=tmpdir,
                        overwrite=overwrite)
        elif records:

            self.sanitise_records()
            if parallel_load:
                self.put_dv_matrices_parallel(helper=helper,
                        tmpdir=tmpdir, overwrite=overwrite)
            else:
                self.put_dv_matrices(helper=helper, tmpdir=tmpdir,
                        overwrite=overwrite)

        self.length = len(self.records)

    def __str__(self):
        s = 'SequenceCollection object:\n'
        s += 'Contains {0} alignments\n'.format(self.length)
        s += 'Source: {0}\n'.format(self.dir)
        return s
