#!/usr/bin/env bash
# Publish the current packaged release to an existing Git repository.

set -euo pipefail

# Require one repository URL or local Git path from the caller.
if [[ $# -ne 1 ]]; then
    printf 'Usage: %s <repository-url>\n' "$0" >&2
    exit 2
fi
repository_url=$1

# Locate this script and read the version declared by the chart program.
script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
version=$(sed -n 's/^__version__ = "\([0-9][0-9.]*\)"$/\1/p' "$script_dir/get_finding_chart.py")
if [[ -z $version ]]; then
    printf 'Could not read the version from get_finding_chart.py.\n' >&2
    exit 1
fi

# Use the matching release folder, or an extracted release folder itself.
release_name="get_finding_chart_$version"
if [[ -d "$script_dir/$release_name" ]]; then
    release_dir="$script_dir/$release_name"
elif [[ $(basename -- "$script_dir") == "$release_name" ]]; then
    release_dir=$script_dir
else
    printf 'Create the %s release folder before updating GitHub.\n' "$release_name" >&2
    exit 1
fi

# Check that the release has its program, documentation, tests, and README image.
for path in get_finding_chart.py README.md CHANGELOG.md LICENSE requirements.txt \
            update_github.sh images/fc_example_g1T.png; do
    if [[ ! -f "$release_dir/$path" ]]; then
        printf 'Release file is missing: %s\n' "$release_dir/$path" >&2
        exit 1
    fi
done
if [[ ! -d "$release_dir/tests" ]]; then
    printf 'Release tests directory is missing: %s\n' "$release_dir/tests" >&2
    exit 1
fi

# Clone the current remote branch into a temporary directory and clean it up on exit.
work_dir=$(mktemp -d)
trap 'rm -rf -- "$work_dir"' EXIT
git clone --quiet -- "$repository_url" "$work_dir/repo"
checkout="$work_dir/repo"

# Copy only release files so local data and old archives never enter the commit.
cp -a -- "$release_dir/get_finding_chart.py" "$release_dir/README.md" \
    "$release_dir/CHANGELOG.md" "$release_dir/LICENSE" \
    "$release_dir/requirements.txt" "$release_dir/update_github.sh" "$checkout/"
mkdir -p -- "$checkout/images" "$checkout/tests"
cp -a -- "$release_dir/images/fc_example_g1T.png" "$checkout/images/"
cp -a -- "$release_dir/tests/." "$checkout/tests/"

# Commit the release files when they changed, then push the cloned branch.
git -C "$checkout" add -- get_finding_chart.py README.md CHANGELOG.md \
    LICENSE requirements.txt update_github.sh images/fc_example_g1T.png tests
if git -C "$checkout" diff --cached --quiet; then
    printf 'GitHub repository already has release v%s.\n' "$version"
    exit 0
fi
git -C "$checkout" commit -m "Release v$version"
git -C "$checkout" push origin HEAD
printf 'Pushed release v%s to %s\n' "$version" "$repository_url"
