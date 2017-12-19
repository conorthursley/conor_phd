#!/bin/bash

# Copyright (C) 2008-2010 Ulf Lorenz
#
# This file is under public domain. You may modify or incorporate it into other
# works without any restrictions.

# Script to be called before running the demos in one block.
#
# Normal users should never require this, it is part of the automatic setup for
# running all demos and inserting the output in the web page.
#
# We do two things here. First, we copy the directory tree to another directory,
# then we enable video and image output in all demos. Copying makes sense insofar,
# as you usually have the sources in your home directory, which is often mounted
# via NFS or something similar. However, the video output is (before postprocessing)
# uncompressed, which can easily result in films of 100 MB per demo. This is sloooow
# if you do this on a network mount.

# Copy the directory tree somewhere else
SOURCE=.
TARGET='/scratch/ulorenz/wp_demos/Demos_4.9.0'


echo "Checking and creating target directory"

if test ! -d $(dirname $TARGET); then
    mkdir $(dirname $TARGET)
fi


echo "Copying files"

cp -pR $SOURCE $TARGET

if test ! $? -eq 0; then
    echo "Copying directory tree failed. Bailing out"
    exit 1
fi

cd $TARGET

# Now turn video/graphics output on. Since this setting is ignored if graphics
# are turned off, we can simply add this to each single initialize.m file.
# Note that we use Windows line ending style. The initialize files are already
# in this format, and this will not change in the nearest future.
echo "Activating movie export"

INIT_FILES=$(find . -name 'qm_init.m')

for f in $INIT_FILES; do
    sed -r 's/^(plots\.density\.type.+)$/\1\nplots.density.export.on = true;\r\nplots.expect.export.on = true;\r/' $f > ttt && mv ttt $f
done
