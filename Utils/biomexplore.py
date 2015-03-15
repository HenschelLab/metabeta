#!/usr/bin/python

from biom.parse import parse_biom_table
import sys

otutable = parse_biom_table(open(sys.argv[1], "U"))
print 'use otutable.iterSamples or otutable.iterObservations, yields triples (abundance-array, name/id, metadata)'
print "OTUs:    %s (%s ...)" % (len(otutable.ObservationIds), ", ".join(otutable.ObservationIds[:3]))
print "samples: %s (%s ...)" % (len(otutable.SampleIds), ", ".join(otutable.SampleIds[:3]))


