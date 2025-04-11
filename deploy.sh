#!/bin/bash

# Get current version
VERSION=$(grep "__version__" wormcat3/__init__.py | cut -d'"' -f2)

IFS='.' read -r -a PARTS <<< "$VERSION"
MAJOR=${PARTS[0]}
MINOR=${PARTS[1]}
PATCH=${PARTS[2]}

# Increment patch
NEW_VERSION="$MAJOR.$MINOR.$((PATCH + 1))"

echo "Bumping version: $VERSION → $NEW_VERSION"

# Update __init__.py
sed -i '' "s/__version__ = .*/__version__ = \"$NEW_VERSION\"/" wormcat3/__init__.py


rm -rf ./dist
rm -rf ./wormcat3.egg-info
python setup.py sdist

# Run twine check and capture the output
CHECK_OUTPUT=$(twine check dist/*)
echo "$CHECK_OUTPUT"

if echo "$CHECK_OUTPUT" | grep -q "PASSED"; then
    echo "Twine check passed, uploading to PyPI..."
    twine upload --repository pypi dist/*
else
    echo "Twine check failed..."
    exit 1
fi
