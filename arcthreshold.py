#! /usr/bin/env python
# encoding: utf-8

__author__ = "asimmons"


class ArcThreshold(object):
    """
    ArcThreshold -
    """

    @staticmethod
    # Should return the same value even if start and end are switched
    def get_string(start, end):
        if str(start) < str(end):
            return str(start) + "_" + str(end)
        else:
            return str(end) + "_" + str(start)
