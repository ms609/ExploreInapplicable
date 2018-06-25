#!/bin/sh

set -e

[ -z "${GITHUB_PAT}" ] && exit 0
[ "${TRAVIS_BRANCH}" != "master" ] && exit 0
git config --global user.email "martins@gmail.com"
git config --global user.name "Martin Smith"
git add --all *
git commit -m"Update compiled files" || true
git push -q
