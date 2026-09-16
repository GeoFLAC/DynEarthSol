#!/bin/sh
# Generate build_revision.hpp: the build-identity facts no predefined macro
# carries, for the build.snapshot block in runtime_info.cxx. Derived every build,
# never hand-written; run from the Makefile, which exports the DES_REV_* values.
#
# Exit status is the interface:
#     0   already correct, nothing written (a no-change make must relink nothing)
#     10  replaced; the caller must drop the object that embeds it
#     *   failed; the caller should stop the build
#
# Rules it must keep:
#   - Degrade, never fail: no git, no tags, no origin, an absent submodule each
#     yield a sentinel. Metadata may not stop a build, let alone a run.
#   - No paths and no credentials out: providers are labels, dirty is counts only
#     (builder=user@host and the origin host are deliberate).
#   - Dimension-neutral: one header serves the 2D and 3D builds.
#   - The "already correct" comparison ignores DES_STATE_UTC, so the clock alone
#     never rewrites the header.
#
# Not in the Makefile: macOS has make 3.81, where .ONESHELL: does not exist and a
# define'd recipe runs one shell per line, so a recipe must be one backslash-joined
# line -- in which a `#` is a make comment that swallows the rest.
# POSIX sh, not python3: gen_feature_deps.py falls back when python3 is missing, so
# it is a soft dependency; this generator has no honest fallback.
# No `set -e`: aborting on the first failed probe is the opposite of rule one.

# Globbing off: git must receive the CODE_PATHS globs unexpanded. Otherwise the
# shell expands them against the current directory and git sees a root-only file
# list, silently dropping nested sources.
set -f

out=$1
[ -n "$out" ] || { echo "usage: sh $0 <output-header>" >&2; exit 2; }

# What the binary is built from -- the only files whose modification can change what
# it computes. Every "is this tree modified" answer here is scoped to them: rev's
# -dirty suffix, the dirty= counts, and the snapshot_diff payload. A README, a .cfg,
# a plot script or a CI workflow cannot change the algorithm. The two build files
# can: they set the per-directory flags the block's resolved compile line does not
# cover (triangle -O1, tetgen -O0, 3x3-C -ffast-math), the source list and the link
# line. benchmarks-cores/Makefile is a test harness, so it is not here. Keep this
# list in step with what the build actually compiles. MMG is the gap the block cannot
# close: its C flags come from its own cmake cache, which no line here records -- the
# mmg= token names the revision in the archive, not what it was compiled with.
CODE_PATHS="*.c *.h *.cxx *.hpp *.cpp *.cu Makefile 3x3-C/Makefile"

# Values make derives with its text functions and exports. Empty is legitimate
# (no warning flags, no GPU_CC); UNSET means the Makefile and this script have
# drifted apart, which is a bug in the caller and must be loud.
for v in DES_REV_DEFINES DES_REV_CXXFLAGS DES_REV_WARNFLAGS DES_REV_LDFLAGS \
         DES_REV_MK_OPTS DES_REV_GPU_CC DES_REV_SNAPSHOT_DIFF \
         DES_REV_KNN_DIR DES_REV_ANN_DIR DES_REV_MMG_DIR \
         DES_REV_PREFIX_BOOST DES_REV_PREFIX_HDF5 DES_REV_PREFIX_MMG \
         DES_REV_PREFIX_OPENMP
do
    eval "set_check=\${$v+set}"
    # shellcheck disable=SC2154  # assigned by the eval above
    [ "$set_check" = set ] || {
        echo "$0: make did not export $v" >&2
        exit 2
    }
done
[ -n "$DES_REV_DEBUG" ] && set -x

# --- dependency providers ----------------------------------------------------
# A prefix path must not reach the binary, so it is reduced to a provider label.
# Not in the Makefile: inside $(shell ...) the ")" closing each `case` pattern
# also closes the expansion.
cls() {
    case "$1" in
        "")                                                   echo toolchain ;;
        /opt/homebrew*|/usr/local/Cellar/*|/usr/local/opt/*)   echo brew ;;
        /opt/local/*)                                          echo macports ;;
        *conda*)                                               echo conda ;;
        "$HOME"*)                                              echo user ;;
        ./mmg/*|mmg/*)                                         echo submodule ;;
        ./external/*|external/*)                               echo vendored ;;
        /opt/nvidia*|*hpc_sdk*)                                echo nvhpc ;;
        *)                                                     echo system ;;
    esac
}
src_boost=$(cls "$DES_REV_PREFIX_BOOST")
src_hdf5=$(cls "$DES_REV_PREFIX_HDF5")
src_mmg=$(cls "$DES_REV_PREFIX_MMG")
src_openmp=$(cls "$DES_REV_PREFIX_OPENMP")

# --- source identity ---------------------------------------------------------
# --tags (version tags are lightweight), --match (a checkpoint tag must not hijack
# the description), --always (a tagless clone still says something). NOT --dirty:
# it judges the whole tree; the suffix comes from the CODE_PATHS diff below.
rev=$(git describe --tags --match "v[0-9]*" --always 2>/dev/null) || rev=unknown
[ -n "$rev" ] || rev=unknown
branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null) || branch=unknown
[ -n "$branch" ] || branch=unknown

# Counts only -- never file names, never content. The sed turns git's prose
# ("2 files changed, 13 insertions(+), 12 deletions(-)") into "2f/+13/-12".
# HEAD-relative so staged counts; additions count here, though the payload omits them.
# An empty shortstat means clean; a FAILED one means nothing was checked, and "clean" is
# a claim rather than the absence of one -- so key on the exit status of this very probe
# (no git, no repository, an unborn HEAD, a locked index all land here).
# shellcheck disable=SC2086  # pathspecs must word-split; set -f keeps them unexpanded
if shortstat=$(git diff --shortstat HEAD -- $CODE_PATHS 2>/dev/null); then
    dirty=$(printf '%s\n' "$shortstat" |
        sed 's/^ *//; s/ files* changed/f/;
             s/, \([0-9]*\) insertions*(+)/\/+\1/;
             s/, \([0-9]*\) deletions*(-)/\/-\1/')
else
    dirty=unknown
fi
[ -n "$dirty" ] && [ "$dirty" != unknown ] && [ "$rev" != unknown ] && rev="$rev-dirty"

# Strip userinfo from both remote spellings (scheme://user@host and the scp form
# user@host:path). A path-only remote names no server, only a directory, so it is
# reduced to a label rather than embedded.
# "none" is a claim about the repository, so keep it for a repository that really has no
# origin: when the probe itself could not run, say unknown, as rev and dirty do.
if url=$(git remote get-url origin 2>/dev/null); then
    origin=$(printf '%s\n' "$url" | sed -e 's|://[^/@]*@|://|' -e 's|^[^/:]*@||')
    case "$origin" in
        /*|.*|file://*) origin=local ;;
    esac
    [ -n "$origin" ] || origin=none
elif git rev-parse --git-dir >/dev/null 2>&1; then
    origin=none      # a repository with no origin remote
else
    origin=unknown   # no git, or not a repository at all
fi

# A submodule directory that is empty is still INSIDE this work tree, so git walks
# up and cheerfully describes the SUPERPROJECT. Require the directory's own .git.
subrev() {
    if [ -n "$1" ] && [ -e "$1/.git" ]; then
        git -C "$1" describe --always --dirty 2>/dev/null || echo unknown
    else
        echo unknown
    fi
}
knnrev=$(subrev "$DES_REV_KNN_DIR")
nfrev=$(subrev "$DES_REV_ANN_DIR")

# The mmg archive is configured once and reused (the Makefile rebuilds it only when
# it is missing), so the checked-out source revision is not necessarily the one in
# the executable. Read the stamp beside the archive; say nothing when there is none.
mmgrev=""
if [ "$src_mmg" = submodule ] && [ -r "$DES_REV_MMG_DIR/build/.mmg-rev" ]; then
    mmgrev="@$(cat "$DES_REV_MMG_DIR/build/.mmg-rev" 2>/dev/null)"
    [ "$mmgrev" = "@" ] && mmgrev=""
fi

# --- build host --------------------------------------------------------------
sutc=$(date -u "+%FT%TZ")
buser=$(id -un 2>/dev/null) || buser=unknown
[ -n "$buser" ] || buser=unknown
bhost=$(hostname -s 2>/dev/null) || bhost=unknown
[ -n "$bhost" ] || bhost=unknown
if command -v sw_vers >/dev/null 2>&1; then
    bos="macos-$(sw_vers -productVersion)"
elif [ -r /etc/os-release ]; then
    # shellcheck disable=SC1091,SC2154  # sourced at run time; ID/VERSION_ID come from it
    bos=$(. /etc/os-release; echo "$ID-$VERSION_ID")
else
    bos=$(uname -s | tr '[:upper:]' '[:lower:]')
fi

# --- emit --------------------------------------------------------------------
# DES_DIRTY carries its own " dirty=" prefix, or is empty on a clean tree, so the
# sentinel's code line needs no conditional. DES_MMG_REV works the same way.
tmp=$out.tmp
{
    printf '#define DES_REVISION "%s"\n'      "$rev"
    printf '#define DES_BRANCH "%s"\n'        "$branch"
    printf '#define DES_DIRTY "%s"\n'         "${dirty:+ dirty=$dirty}"
    printf '#define DES_ORIGIN "%s"\n'        "$origin"
    printf '#define DES_STATE_UTC "%s"\n'     "$sutc"
    printf '#define DES_BUILD_OS "%s"\n'      "$bos"
    printf '#define DES_BUILDER "%s"\n'       "$buser@$bhost"
    printf '#define DES_MAKE_OPTS "%s"\n'     "$DES_REV_MK_OPTS"
    printf '#define DES_GPU_CC "%s"\n'        "$DES_REV_GPU_CC"
    printf '#define DES_KNNBVH_REV "%s"\n'    "$knnrev"
    printf '#define DES_NANOFLANN_REV "%s"\n' "$nfrev"
    printf '#define DES_MMG_REV "%s"\n'       "$mmgrev"
    printf '#define DES_SRC_BOOST "%s"\n'     "$src_boost"
    printf '#define DES_SRC_HDF5 "%s"\n'      "$src_hdf5"
    printf '#define DES_SRC_MMG "%s"\n'       "$src_mmg"
    printf '#define DES_SRC_OPENMP "%s"\n'    "$src_openmp"
    printf '#define DES_DEFINES "%s"\n'       "$DES_REV_DEFINES"
    printf '#define DES_CXXFLAGS "%s"\n'      "$DES_REV_CXXFLAGS"
    printf '#define DES_WARNFLAGS "%s"\n'     "$DES_REV_WARNFLAGS"
    printf '#define DES_LDFLAGS "%s"\n'       "$DES_REV_LDFLAGS"
} > "$tmp" || { echo "$0: cannot write $tmp" >&2; exit 1; }

# --- the uncommitted code changes, only when asked for ----------------------
# Off by default for privacy. Compiled sources only (CODE_PATHS), tracked only
# (--diff-filter=a drops additions, so a file absent from HEAD contributes nothing,
# not even its name), HEAD-relative; whatever is left out is counted below.
if [ "$DES_REV_SNAPSHOT_DIFF" = 1 ]; then
    body=$out.body
    {
        echo '    '; echo '==== Summary of the code ===='; echo '    '
        if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
            git show -s 2>&1; echo '    '
            # shellcheck disable=SC2086  # must word-split into pathspecs (set -f above
            # keeps the shell from expanding them; git needs the globs verbatim)
            git diff --stat HEAD --diff-filter=a -- $CODE_PATHS 2>&1; echo '    '
            # Name the excluded counts, never the files: a path is a path. Both
            # excluded classes are disclosed -- non-code files, and added sources.
            nall=$(git diff HEAD --name-only --diff-filter=a 2>/dev/null | wc -l | tr -d ' ')
            # shellcheck disable=SC2086  # pathspecs must word-split; set -f keeps them unexpanded
            ncode=$(git diff HEAD --name-only --diff-filter=a -- $CODE_PATHS 2>/dev/null | wc -l | tr -d ' ')
            # shellcheck disable=SC2086
            nadd=$(git diff HEAD --name-only --diff-filter=A -- $CODE_PATHS 2>/dev/null | wc -l | tr -d ' ')
            if [ "${nall:-0}" -gt "${ncode:-0}" ]; then
                echo "   $((nall - ncode)) non-code file(s) also changed; not embedded (only $CODE_PATHS are)"
            fi
            if [ "${nadd:-0}" -gt 0 ]; then
                echo "   $nadd added source file(s) counted in dirty= but not embedded (not in HEAD)"
            fi
            if [ "${nall:-0}" -gt "${ncode:-0}" ] || [ "${nadd:-0}" -gt 0 ]; then
                echo '    '
            fi
            echo '== Code modification (not checked-in) =='; echo ' '
            # shellcheck disable=SC2086  # pathspecs must word-split; set -f keeps them unexpanded
            if [ -n "$(git diff HEAD --name-only --diff-filter=a -- $CODE_PATHS 2>/dev/null)" ]; then
                # shellcheck disable=SC2086
                git diff HEAD --diff-filter=a -- $CODE_PATHS 2>&1
            else
                # An embedded-but-empty payload must not read like a missing one.
                echo '   (none: no tracked source differs from HEAD)'
            fi
            echo ' '
            # Skip the section rather than embed git's error text when the branch
            # has no upstream and there is no origin/HEAD to fall back on.
            up=$(git rev-parse --verify --quiet '@{upstream}' || git rev-parse --verify --quiet origin/HEAD)
            if [ -n "$up" ]; then
                echo '== Commits not in upstream =='; echo '    '
                git log --oneline "$up".. 2>&1; echo '    '
            fi
        else
            echo 'git or the repository is unavailable: no history or diff recorded'
            echo '    '
        fi
    } > "$body"

    # 1 MB, far above any development diff, measured before the indent below.
    if [ "$(wc -c < "$body")" -gt 1000000 ]; then
        { head -c 1000000 "$body"; printf '\n[truncated at 1 MB]\n'; } > "$body.2"
        mv "$body.2" "$body"
    fi

    # A raw string needs a delimiter the payload does not contain -- and a diff of
    # THIS FILE contains any fixed choice, including this comment. Search for one.
    d=DESDIFF
    n=0
    while grep -q ")$d\"" "$body"; do
        n=$((n + 1))
        d=DESDIFF$n
    done

    {
        echo '// code changes embedded by snapshot_diff=1'
        echo '#define DES_HAS_CODE_DIFF 1'
        echo '__attribute__((used)) const char snapshot_code_diff[] ='
        printf '    "\\nbuild.code-changes.begin :"\n'
        printf 'R"%s(\n' "$d"
        sed '/^== Code modification (not checked-in) ==/,/^== Commits not in upstream ==/{/^==/!s/^/   /;}' "$body"
        printf ')%s"\n' "$d"
        printf '    "build.code-changes.end   :";\n'
    } >> "$tmp"
else
    echo '// snapshot_diff=0: no code changes embedded (make snapshot_diff=1 to embed).' >> "$tmp"
fi
rm -f "$out.body"

# --- replace only on a real change ------------------------------------------
# state_utc is excluded from the comparison: it is the only field the clock moves,
# and letting it force a rewrite would relink the world on every make.
grep -v '^#define DES_STATE_UTC ' "$tmp" > "$tmp.cmp1"
grep -v '^#define DES_STATE_UTC ' "$out" > "$tmp.cmp2" 2>/dev/null || :
if cmp -s "$tmp.cmp1" "$tmp.cmp2"; then
    rm -f "$tmp" "$tmp.cmp1" "$tmp.cmp2"
    exit 0
fi
mv "$tmp" "$out" || { rm -f "$tmp" "$tmp.cmp1" "$tmp.cmp2"; exit 1; }
rm -f "$tmp.cmp1" "$tmp.cmp2"
exit 10
