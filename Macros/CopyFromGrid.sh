#!/bin/bash
# Load O2 env before running this script
TARGET_DIR=$1
LOCAL_DIR=$2
#RUN=$3
#PASS=$4

# NEW_DIR=$PASS
# mkdir -p $LOCAL_DIR/$NEW_DIR

cd $LOCAL_DIR/
# mkdir -p $RUN_DEBUG
# cd $RUN_DEBUG
echo "Start processing $TARGET_DIR..."

# COPY_DIR=$TARGET_DIR/
# echo "Start copying from $COPY_DIR ..."
# CTF_LIST=`alien_ls $TARGET_DIR`
CTF_LIST=`alien_ls $TARGET_DIR | grep o2_ctf*`
echo "Found jobs:"
echo $CTF_LIST
for CTF in $CTF_LIST
do
	mkdir -p $CTF
	echo "Processing $CTF"

	for ITEM in mchtracks.root muontracks.root workflowconfig.log
	do
		echo "Copying $ITEM ..."
		alien_cp alien:$TARGET_DIR/$CTF/$ITEM file:$LOCAL_DIR/$CTF/$ITEM
		echo "$ITEM copied."
	done

done

cd $LOCAL_DIR/
# hadd -f mchtracks.root o2_ctf*/mchtracks.root
# hadd -f muontracks.root o2_ctf*/muontracks.root
# hadd -f mid-reco.root o2_ctf*/mid-reco.root