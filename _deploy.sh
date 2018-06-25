#!/bin/sh

set -e

[ -z "${GITHUB_PAT}" ] && exit 0
[ "${TRAVIS_BRANCH}" != "master" ] && exit 0
git config --global user.email "martins@gmail.com"
git config --global user.name "Martin Smith"
git clone -b origin https://${GITHUB_PAT}@github.com/${TRAVIS_REPO_SLUG}.git git-files
cd git-files
cp -r ../*.html ./
git add --all *
git commit -m"[TRAVIS:] Update compiled html files" || true
git push -q origin
