#!/usr/bin/sh
cat pyproject.toml | grep "version = \"[^\"]*\"" | grep -oP "\"[^\"]*\"" | grep -oP "[^\"]*"  