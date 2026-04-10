#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <source.txt> <target string> <anchor string>" >&2
    exit 1
fi

src=$1
target=$2
anchor=$3

if [ ! -f "$src" ]; then
    echo "Source file not found: $src" >&2
    exit 1
fi

tmp=$(mktemp) || exit 1
trap 'rm -f "$tmp"' EXIT

LC_ALL=en_GB.UTF-8 awk -v target="$target" -v anchor="$anchor" '
BEGIN {
    target_found = 0
    anchor_found = 0
}
{
    lines[++n] = $0

    if (!target_found && index($0, target) > 0) {
        target_line = n
        target_found = 1
    }
}
END {
    if (!target_found) {
        print "Target string not found in any line." > "/dev/stderr"
        exit 1
    }

    m = 0
    for (i = 1; i <= n; i++) {
        if (i != target_line) {
            filtered[++m] = lines[i]
        }
    }

    for (i = 1; i <= m; i++) {
        if (!anchor_found && index(filtered[i], anchor) > 0) {
            anchor_line = i
            anchor_found = 1
        }
    }

    if (!anchor_found) {
        print "Anchor string not found in any line." > "/dev/stderr"
        exit 1
    }

    insert_line = anchor_line - 3
    if (insert_line < 1) {
        insert_line = 1
    }

    for (i = 1; i < insert_line; i++) {
        print filtered[i]
    }

    print target

    for (i = insert_line; i <= m; i++) {
        print filtered[i]
    }
}
' "$src" > "$tmp"

mv "$tmp" "$src"
trap - EXIT

echo "Done."