#!/usr/bin/python
# -*- coding: utf-8 -*-


class Result(object):

    def update_score(self):
        self.score = sum([rec.tree.score for rec in self.concats])

    def __init__(self, clusters={}):
        self.score = 0
        self.length = 0
        self.concats = []
        self.members = []
        self.names = []
        if clusters:
            self.length = len(clusters)
            for key in sorted(clusters):
                self.concats.append(clusters[key]['concatenation'])
                self.members.append(clusters[key]['members'])
                self.names.append([rec.name for rec in
                                  clusters[key]['members']])
            self.update_score()

    def retrieve_names(self, index):
        return self.names[index]

    def retrieve_concat(self, index):
        return self.concats[index]

    def retrieve_members(self, index):
        return self.members[index]
