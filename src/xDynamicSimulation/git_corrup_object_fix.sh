#!/bin/bash 

#specify SIRF git repository and the forked repository
remote_url_SIRF="https://github.com/CCPPETMR/SIRF.git"
remote_url_FORK="https://github.com/johannesmayer/SIRF.git"

# specify branches
dev_branch="devJohannes"
curr_branch="noncartMrAcquModel"

# generate a hard backup for SIRF
todays_date=`date +%Y-%m-%d`
backup_filepath="/home/sirfuser/Code/Backup/SIRF_BACKUP_$todays_date"

echo "Generating $backup_filepath for backing up SIRF"

mkdir -p $backup_filepath

echo "Copying SIRF into the backup filepath"
cp -r $SIRF_PATH $backup_filepath


path_to_sirf="/home/sirfuser/devel/buildVM/sources/"
cd $path_to_sirf

rm -rf $SIRF_PATH


if [ ! -d "$SIRF_PATH" ];then
	echo "Cloning SIRF repository"
	git clone $remote_url_SIRF
fi

cd $SIRF_PATH

git remote -v

echo "Temporarily setting the origin to the fork repository"
git remote set-url origin $remote_url_FORK

git remote -v

echo "Pulling branches from fork"

git pull origin master
git fetch
git checkout -b $dev_branch
git pull origin $dev_branch

git fetch 
git checkout -b $curr_branch
git pull origin $curr_branch

echo "Setting back the pull to SIRF and the push to FORK"

git remote set-url origin $remote_url_SIRF
git remote set-url --push origin $remote_url_FORK

git remote -v


echo "Copying current file state from the backup back to SIRF"
cp -r "$backup_filepath/SIRF/src/xDynamicSimulation" "$SIRF_PATH/src"



