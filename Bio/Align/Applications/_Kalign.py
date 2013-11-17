# Copyright 2013 by xyz.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program Kalign.
"""

__docformat__ = "epytext en" #Don't just use plain text in epydoc API pages!

import os
from Bio.Application import _Option, _Switch, AbstractCommandline

class KalignCommandline(AbstractCommandline):
    """Command line wrapper for Kalign.

    http://msa.sbc.su.se/cgi-bin/msa.cgi -> Kalign

    Example:

    >>> from Bio.Align.Applications import KalignCommandline
    >>> kalign_cline = KalignCommandline(input="kalign_in.fasta", \
                                         output="kalign_out.fasta")
    >>> print kalign_cline
    kalign -input=kalign_in.fasta -output=kalign_out.fasta

    You would typically run the command line with kalign_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Last checked against version: 2.04
    """
    def __init__(self, cmd="kalign", **kwargs):
        self.parameters = \
            [
            _Option(["-input", "input"],
                    "Name of the input file.",
                    filename=True, is_required=True),
            _Option(["-output", "output"],
                    "Name of the output file.",
                    filename=True, is_required=True),
            _Option(["-gpo", "gpo"],
                    "Gap open penalty.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-gpe", "gpe"],
                    "Gap extension penalty.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-tgpe", "tgpe"],
                    "Terminal gap penalties.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-bonus", "bonus"],
                    "A constant added to the substitution matrix.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-sort", "sort"],
                    "The order in which the sequences appear in the output"
                    " alignment.",
                    checker_function=lambda x: x in ["input", "tree", "gaps"]),
            _Option(["-feature", "feature"],
                    "Selects feature mode and specifies which features are"
                    " to be used. Valid values: all, maxplp, STRUCT, PFAM-A.",
                    checker_function=lambda x: x in ["all", "maxplp",
                                                     "STRUCT", "PFAM-A"]),
            _Option(["-same_feature_score", "same_feature_score"],
                    "Score for aligning same features.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-diff_feature_score", "diff_feature_score"],
                    "Penalty for aligning different features.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-distance", "distance"],
                    "Distance method",
                    checker_function=lambda x: x in ["wu", "pair"]),
            _Option(["-tree", "tree"],
                    "Guide tree method",
                    checker_function=lambda x: x in ["nj", "upgma"]),
            _Option(["-zcutoff", "zcutoff"],
                    "Parameter used in the wu-manber based distance"
                    " calculation.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-gap_inc", "gap_inc"],
                    "Increases gap penalties depending on the number"
                    " of existing gaps.",
                    checker_function=lambda x: isinstance(x, int) or \
                                               isinstance(x, float)),
            _Option(["-format", "format"],
                    "Output format.",
                    checker_function=lambda x: x in ["fasta", "msf", "aln",
                                                     "clu", "macsim"]),
            _Switch(["-quiet", "quiet"],
                    "Print nothing to STDERR. Read nothing from STDIN."),
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

def _test():
    """Run the module's doctests (PRIVATE)."""
    print "Runing Kalign doctests..."
    import doctest
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    _test()
