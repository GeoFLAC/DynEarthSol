#!/bin/sh
# Generate build_revision.hpp, the build-identity macros for the build.snapshot block in
# runtime_info.cxx; the Makefile runs it and exports the DES_REV_* values. Exit 0: already
# correct, nothing written; 10: replaced, so drop the object embedding it; else: failed.
#
# Rules: degrade, never fail (no git, tag, origin or submodule yields a sentinel); no paths
# or credentials out (providers are labels, dirty is counts; builder= and the origin host
# are deliberate); one header for 2D and 3D; the clock alone never rewrites the header.
#
# A script, not a recipe: macOS's make 3.81 runs a recipe one shell per line. POSIX sh,
# not python3, which is optional in this build. No `set -e`: see the first rule.

# git must get the CODE_PATHS globs unexpanded, or nested sources silently drop out.
set -f

out=$1
[ -n "$out" ] || { echo "usage: sh $0 <output-header>" >&2; exit 2; }

# The files that can change what the binary computes, sources and the two build files
# (they set per-directory flags like triangle -O1 that the compile line misses); every
# "modified" answer here is scoped to them. MMG's C flags live in its cmake cache, so its
# mmg= token names the revision only.
CODE_PATHS="*.c *.h *.cxx *.hpp *.cpp *.cu Makefile 3x3-C/Makefile"

# Empty is legitimate (no GPU_CC); UNSET means the Makefile and this script drifted apart.
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
# A prefix path must not reach the binary: it becomes a provider label.
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
# --tags: version tags are lightweight; --match: no checkpoint tag hijacks it; --always:
# a tagless clone still answers. No --dirty, which judges the whole tree.
rev=$(git describe --tags --match "v[0-9]*" --always 2>/dev/null) || rev=unknown
[ -n "$rev" ] || rev=unknown
branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null) || branch=unknown
[ -n "$branch" ] || branch=unknown

# Counts only, never names: "2 files changed, 13 insertions(+), 12 deletions(-)" becomes
# "2f/+13/-12", HEAD-relative. A failed probe checked nothing, so it reads unknown, not clean.
# shellcheck disable=SC2086
if shortstat=$(git diff --shortstat HEAD -- $CODE_PATHS 2>/dev/null); then
    dirty=$(printf '%s\n' "$shortstat" |
        sed 's/^ *//; s/ files* changed/f/;
             s/, \([0-9]*\) insertions*(+)/\/+\1/;
             s/, \([0-9]*\) deletions*(-)/\/-\1/')
else
    dirty=unknown
fi
[ -n "$dirty" ] && [ "$dirty" != unknown ] && [ "$rev" != unknown ] && rev="$rev-dirty"

# Strip userinfo from both remote spellings; a path-only remote becomes "local". "none" is
# a repository without origin; a probe that could not run says unknown.
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

# In an empty submodule directory git would describe the superproject: require its .git.
subrev() {
    if [ -n "$1" ] && [ -e "$1/.git" ]; then
        git -C "$1" describe --always --dirty 2>/dev/null || echo unknown
    else
        echo unknown
    fi
}
knnrev=$(subrev "$DES_REV_KNN_DIR")
nfrev=$(subrev "$DES_REV_ANN_DIR")

# The stamp beside the mmg archive (the Makefile writes it) names what links in.
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
# DES_DIRTY and DES_MMG_REV carry their own prefix, or are empty: no conditional needed.
tmp=$out.tmp.$$   # per process: two builds in one tree may run this at once
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
# Off by default for privacy: tracked CODE_PATHS only (--diff-filter=a drops added files,
# names included), HEAD-relative; what is left out is counted.
if [ "$DES_REV_SNAPSHOT_DIFF" = 1 ]; then
    body=$tmp.body
    {
        echo '    '; echo '==== Summary of the code ===='; echo '    '
        if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
            git show -s 2>&1; echo '    '
            # shellcheck disable=SC2086
            git diff --stat HEAD --diff-filter=a -- $CODE_PATHS 2>&1; echo '    '
            # Count what is left out, never name it.
            nall=$(git diff HEAD --name-only --diff-filter=a 2>/dev/null | wc -l | tr -d ' ')
            # shellcheck disable=SC2086
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
            # shellcheck disable=SC2086
            if [ -n "$(git diff HEAD --name-only --diff-filter=a -- $CODE_PATHS 2>/dev/null)" ]; then
                # shellcheck disable=SC2086
                git diff HEAD --diff-filter=a -- $CODE_PATHS 2>&1
            else
                # An embedded-but-empty payload must not read like a missing one.
                echo '   (none: no tracked source differs from HEAD)'
            fi
            echo ' '
            # Without an upstream or origin/HEAD, skip the section rather than embed an error.
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
rm -f "$tmp.body"

# --- replace only on a real change ------------------------------------------
# Compared without DES_STATE_UTC, or every make would relink everything.
grep -v '^#define DES_STATE_UTC ' "$tmp" > "$tmp.cmp1"
grep -v '^#define DES_STATE_UTC ' "$out" > "$tmp.cmp2" 2>/dev/null || :
if cmp -s "$tmp.cmp1" "$tmp.cmp2"; then
    rm -f "$tmp" "$tmp.cmp1" "$tmp.cmp2"
    exit 0
fi
mv "$tmp" "$out" || { rm -f "$tmp" "$tmp.cmp1" "$tmp.cmp2"; exit 1; }
rm -f "$tmp.cmp1" "$tmp.cmp2"
exit 10
