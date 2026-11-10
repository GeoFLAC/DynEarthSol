#!/bin/bash
# Deletes stale versions of the ghcr.io layer cache from docker.yml:
#   - untagged versions, replaced by a later publish;
#   - versions tagged only -pr<N> of closed pull requests.
# A version that also has master's tag is kept.
#
# DRY_RUN=1 only lists. Needs GH_TOKEN (packages: write, pull-requests: read),
# GH_REPO, and the repository's Admin role on the package, else 403.
set -eo pipefail

API=/orgs/GeoFLAC/packages/container/dynearthsol%2Fdocker-cache/versions

# Per version: "<id> <pr-or-dash per tag>", e.g. "7 - 119".
versions=$(gh api --paginate "$API?per_page=100" --jq '.[] |
    [(.id | tostring)] + [.metadata.container.tags[] |
        (capture("-pr(?<n>[0-9]+)$").n // "-")] | join(" ")')

open=" $(gh pr list --state open --limit 1000 --json number -q '.[].number' | tr '\n' ' ') "

ids=
while read -r id prs; do
    [ -n "$id" ] || continue
    keep=
    for n in $prs; do
        [[ $n = - || $open == *" $n "* ]] && { keep=1; break; }
    done
    [ -z "$keep" ] && ids="$ids $id"
done <<< "$versions"

echo "$(echo "$ids" | wc -w | tr -d ' ') version(s) to delete"
for id in $ids; do
    if [ "$DRY_RUN" = 1 ]; then
        echo "would delete $id"
    else
        # 404 means already deleted; any other error stops.
        if out=$(gh api -X DELETE "$API/$id" 2>&1); then
            echo "deleted $id"
        elif [[ $out == *"HTTP 404"* ]]; then
            echo "already gone $id"
        else
            echo "$out"; exit 1
        fi
    fi
done
