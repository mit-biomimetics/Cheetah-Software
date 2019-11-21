#!/bin/bash
set -ev


echo "Deploying docs to website osqp.org..."

# Update variables from install

# Anaconda
export PATH=${DEPS_DIR}/miniconda/bin:$PATH
hash -r
source activate testenv


if [[ $TRAVIS_OS_NAME == "linux" ]]; then

	# Install dependencies
	# Add PPA to make doxygen installation work
	sudo add-apt-repository ppa:libreoffice/ppa -y
	sudo apt-get update -q -y
	sudo apt-get install -y doxygen
	conda install -y sphinx
	conda install -y -c conda-forge sphinx_rtd_theme breathe

	# Enter in docs folder
	cd ${TRAVIS_BUILD_DIR}/docs/

	# Create docs
	make html

	# Clone website repository
	cd $HOME/
	git clone "https://${OSQP_DOCS_DEPLOY_GH_TOKEN}@github.com/oxfordcontrol/osqp.org.git"
	cd osqp.org/src/docs/

	# Update git config to push
	git config user.name "OSQP Docs Builder"
	git config user.email "$OSQP_DOCS_DEPLOY_EMAIL"

	# Copy in the docs HTML
	cp -R ${TRAVIS_BUILD_DIR}/docs/_build/html/* ./

	# Add and commit changes.
	git add -A .
	git commit --allow-empty -m "Generated docs for commit $TRAVIS_COMMIT"
	# -q is very important, otherwise you leak your GH_TOKEN
	git push -q origin master

fi

