#!/bin/bash

# Copyright (C) 2008-2010 Ulf Lorenz
#
# This file is under public domain. You may modify or incorporate it into other
# works without any restrictions.

# Script to be called after running the demos in one block.
#
# Normal users should never require this, it is part of the automatic setup for
# running all demos and inserting the output in the web page.
#
# We have now copied the files to another path, setup the initialize.m files,
# and run all the demos. Now we have to clean up. There are two things to do.

# the directory to process
SOURCE=.

cd $SOURCE

# 1. The movies are uncompressed right now. We have to convert them to something
#    else (we use mpeg). See user reference, I just copy the entry there.
SURFACE_FILMS=$(find . -name 'surface.avi' | sort -g)
CONTOUR_FILMS=$(find . -name 'contour.avi' | sort -g)
REDUCED_FILMS=$(find . -name 'reduced.avi' | sort -g)
CURVE_FILMS=$(find . -name 'curve.avi' | sort -g)
POLAR_FILMS=$(find . -name 'polar.avi' | sort -g)


for f in $SURFACE_FILMS $CONTOUR_FILMS $REDUCED_FILMS $CURVE_FILMS $POLAR_FILMS; do
    echo "Processing $f"
    dir=$(dirname $f)
    file=$(basename $f .avi)

    pushd $dir &>/dev/null
	ffmpeg -i ${file}.avi -vcodec mpeg1video -r 25 -b:v 1500k -filter:v "setpts=3.0*PTS" -y ${file}.mpg

	# -----------------------------------------------------------------------
	# Alternative encoding:
	#
	# Re-encoding the film in mpeg format has the disadvantage that mpeg was
	# meant for "real-life" data. That is, it works pretty poorly for the plots
	# we have here, with very steep color gradients (from background to drawn
	# line). As a result, you will often see block artefacts, especially if
	# curves change a lot. Also, the gray background gets a definite pinkish
	# touch.
    #
    # See the FAQ on the Wiki for how to use other formats, such as qtrle, for example.
	# -----------------------------------------------------------------------

    if test ! $? -eq 0; then
        echo "Error converting film $f"
        exit 1
    fi

    popd &>/dev/null
done

# 2. Remove all files that are not needed. This includes the uncompressed films,
#    the subversion directories, and all the automatic runscripts.
find . -name '*.avi' -delete
find . -path '*.svn*' -delete

# The remaining scripts you can easily delete yourself.
