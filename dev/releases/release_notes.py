#!/usr/bin/env python3
#############################################################################
# Usage:
#     ./release_notes.py [VERSION]
#
# For example
#     ./release_notes.py 4.13.1
#
# This assumes that the tags named v4.13.1, 4.13dev (?) and v4.13.0 (???) already exists.
#
# A version ending in .0 is consider MAJOR, any other MINOR
# Don't use this with versions like 4.13.0-beta1

import json
import os
import re
import copy
import subprocess
import sys
from datetime import datetime
from typing import Any, Dict, List
import tomli

ownpath = os.path.abspath(sys.argv[0])
dirpath = os.path.dirname(ownpath)
repopath = os.path.dirname(os.path.dirname(os.path.dirname(ownpath)))
NEWFILE = f"{dirpath}/new.md"
FINALFILE = f"{repopath}/CHANGELOG.md"

# read config file
with open('config.toml', 'rb') as conffile:
    conf = tomli.load(conffile)

REPONAME = conf['reponame']
PROJECTNAME = REPONAME.split('/')[-1]
ENABLE_TWOLEVEL = conf["enabletwolevel"]
REQUIRE_TWOLEVEL = False
if ENABLE_TWOLEVEL:
    REQUIRE_TWOLEVEL = conf['requiretwolevel']

# the following loads a dict of {LABEL: DESCRIPTION}; the first entry is the name of a GitHub label
# (be careful to match them precisely), the second is a headline for a section the release notes;
# any PR with the given label is put into the corresponding section; each PR is put into only one
# section, the first one from this list it fits in.
# See also <https://github.com/gap-system/gap/issues/4257>.

TOPICS = conf['topics']
PRTYPES = {}
if ENABLE_TWOLEVEL:
    PRTYPES = conf['prtypes']


def usage(name: str) -> None:
    print(f"Usage: `{name} [NEWVERSION]`")
    sys.exit(1)


def is_existing_tag(tag: str) -> bool:
    print(tag)
    res = subprocess.run(
        [
            "gh",
            "release",
            "list",
            "--json=name",
            "-q",
            f""".[] | select(.name | contains("{tag.strip()}"))"""
        ],
        shell=False,
        check=False, # this subprocess is allowed to fail
        capture_output=True
    )
    return res.stdout.decode() != ""


def find_previous_version(version: str) -> str:
    major, minor, patchlevel = map(int, version.split("."))
    if patchlevel != 0:
        patchlevel -= 1
        return f"{major}.{minor}.{patchlevel}"
    minor -= 1
    patchlevel = 0
    while True:
        v = f"{major}.{minor}.{patchlevel}"
        if not is_existing_tag("v" + v):
            break
        patchlevel += 1
    if patchlevel == 0:
        error("could not determine previous version")
    patchlevel -= 1
    return f"{major}.{minor}.{patchlevel}"

def notice(s):
    print(s)

def error(s):
    print(s)
    sys.exit(1)

def warning(s):
    print('===================================================')
    print(s)
    print('===================================================')



def get_tag_date(tag: str) -> str:
    if is_existing_tag(tag):
        res = subprocess.run(
            [
                "gh",
                "release",
                "view",
                f"{tag}",
                "--json=createdAt"
            ],
            shell=False,
            check=True,
            capture_output=True
        )
        res = json.loads(res.stdout.decode())
    else:
        error("tag does not exist!")
    return res['createdAt'][0:10]


def get_pr_list(date: str, extra: str) -> List[Dict[str, Any]]:
    query = (
        f'merged:>={date} -label:"release notes: not needed" -label:"release notes: added"'
        f'base:master {extra}'
    )
    print("query: ", query)
    res = subprocess.run(
        [
            "gh",
            "pr",
            "list",
            "--search",
            query,
            "--json",
            "number,title,closedAt,labels,mergedAt,body",
            "--limit",
            "200",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    json_list = json.loads(res.stdout.strip())
    json_list = sorted(json_list, key=lambda d: d['number']) # sort by ascending PR number
    return json_list


def pr_to_md(pr: Dict[str, Any]) -> str:
    """Returns markdown string for the PR entry"""
    k = pr["number"]
    if has_label(pr, 'release notes: use body'):
        mdstring = re.sub(
            r'^- ', f"- [#{k}](https://github.com/{REPONAME}/pull/{k}) ",
            pr["body"]
        )
    else:
        title = pr["title"]
        mdstring = f"- [#{k}](https://github.com/{REPONAME}/pull/{k}) {title}\n"
    return mdstring

def body_to_release_notes(pr):
    body = pr['body']
    index1 = body.lower().find("## release notes")
    if index1 == -1:
        ## not found
        ## complain and return fallback
        print(f"Release notes section not found in PR number {pr['number']}!!")
        return body
    index2 = body.find('\n', index1) + 1 # the first line after the release notes line
    bodylines = body[index2:].splitlines()
    mdstring = ""
    for line in bodylines:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('- '):
            mdstring = f"{mdstring}\n{line}"
        else:
            break
    if not mdstring:
        warning(f"Empty release notes section for PR #{pr['number']} !")
    return mdstring


def has_label(pr: Dict[str, Any], label: str) -> bool:
    return any(x["name"] == label for x in pr["labels"])


def changes_overview(
    prs: List[Dict[str, Any]], new_version: str
) -> None:
    """Writes files with information for release notes."""

    date = datetime.now().strftime("%Y-%m-%d")
    release_url = f"https://github.com/{REPONAME}/releases/tag/v{new_version}"

    # Could also introduce some consistency checks here for wrong combinations of labels
    notice("Writing release notes into file " + NEWFILE)
    with open(NEWFILE, "w", encoding="utf-8") as relnotes_file:
        prs_with_use_title = [
            pr for pr in prs if
            has_label(pr, "release notes: use title") or
            has_label(pr, "release notes: use body")
        ]
        # Write out all PRs with 'use title'
        relnotes_file.write(
            f"""# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [{new_version}]({release_url}) - {date}

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

"""
        )
        totalPRs = len(prs)
        print(f"Total number of PRs: {totalPRs}")
        countedPRs = 0
        for priorityobject in TOPICS:
            matches = [
                pr for pr in prs_with_use_title if has_label(pr, priorityobject)
            ]
            original_length = len(matches)
            print("PRs with label '" + priorityobject + "': ", len(matches))
            countedPRs = countedPRs + len(matches)
            if len(matches) == 0:
                continue
            relnotes_file.write("### " + TOPICS[priorityobject] + "\n\n")
            if priorityobject == "breaking":
                relnotes_file.write("> !These changes break compatibility from previous versions!")
                relnotes_file.write("\n\n")
            if priorityobject in ['release notes: highlight', 'breaking']:
                itervar = TOPICS
            else:
                itervar = PRTYPES
            for typeobject in itervar:
                if typeobject == priorityobject:
                    continue
                matches_type = [
                    pr for pr in matches if has_label(pr, typeobject)
                ]
                print(
                    f"PRs with label '{priorityobject}' and type '{typeobject}': "
                    f"{len(matches_type)}"
                )
                if len(matches_type) == 0:
                    continue
                relnotes_file.write(f"#### {itervar[typeobject]}\n\n")
                for pr in matches_type:
                    relnotes_file.write(pr_to_md(pr))
                    prs_with_use_title.remove(pr)
                    matches.remove(pr)
                relnotes_file.write('\n')
            # Items without a type label
            if len(matches) > 0:
                if len(matches) != original_length:
                    relnotes_file.write("#### Miscellaneous changes\n\n")
                for pr in matches:
                    relnotes_file.write(pr_to_md(pr))
                    prs_with_use_title.remove(pr)
                relnotes_file.write('\n')

        print(f"Remaining PRs: {totalPRs - countedPRs}")
        # The remaining PRs have no "kind" or "topic" label from the priority list
        # (may have other "kind" or "topic" label outside the priority list).
        # Check their list in the release notes, and adjust labels if appropriate.
        if len(prs_with_use_title) > 0 and ENABLE_TWOLEVEL:
            relnotes_file.write("### Other changes\n\n")
            for typeobject in PRTYPES:
                matches_type = [
                    pr for pr in prs_with_use_title if has_label(pr, typeobject)
                ]
                len(matches_type)
                print("PRs with type '" + typeobject + "': ", len(matches_type))
                if len(matches_type) == 0:
                    continue
                relnotes_file.write("#### " + PRTYPES[typeobject] + "\n\n")

                for pr in matches_type:
                    relnotes_file.write(pr_to_md(pr))
                    prs_with_use_title.remove(pr)
                relnotes_file.write("\n")

        # Report PRs that have to be updated before inclusion into release notes.
        prs_to_be_added = [pr for pr in prs if has_label(pr, "release notes: to be added")]
        if len(prs_to_be_added) > 0:
            relnotes_file.write("### **TODO** release notes: to be added" + "\n\n")
            relnotes_file.write(
                "If there are any PRs listed below, check their title and labels.\n"
            )
            relnotes_file.write(
                'When done, change their label to "release notes: use title".\n\n'
            )
            for pr in prs_to_be_added:
                relnotes_file.write(pr_to_md(pr))
            relnotes_file.write("\n")
        if len(prs_with_use_title) > 0:
            relnotes_file.write(
                "### **TODO** insufficient labels for automatic classification\n\n"
                "The following PRs have neither a topic label assigned to them, nor a PR type. \n"
                "**Manual intervention required.**\n\n")
            for pr in prs_with_use_title:
                relnotes_file.write(pr_to_md(pr))
                relnotes_file.write('\n')
            relnotes_file.write('\n')

        # remove PRs already handled earlier
        prs = [pr for pr in prs if not has_label(pr, "release notes: to be added")]
        prs = [pr for pr in prs if not has_label(pr, "release notes: added")]
        prs = [pr for pr in prs if not has_label(pr, "release notes: use title")]
        prs = [pr for pr in prs if not has_label(pr, "release notes: use body")]

        # Report PRs that have neither "to be added" nor "added" or "use title" label
        if len(prs) > 0:
            relnotes_file.write("### **TODO** Uncategorized PR" + "\n\n")
            relnotes_file.write(
                "If there are any PRs listed below, either apply the same steps\n"
            )
            relnotes_file.write(
                'as above, or change their label to "release notes: not needed".\n\n'
            )
            for pr in prs:
                relnotes_file.write(pr_to_md(pr))
            relnotes_file.write('\n')

        # now read back the rest of changelog.md into newfile
        with open(FINALFILE, 'r', encoding='ascii') as oldchangelog:
            oldchangelog.seek(262)
            for line in oldchangelog.readlines():
                relnotes_file.write(line)
        # finally copy over this new file to changelog.md
        os.rename(NEWFILE, FINALFILE)

def split_pr_into_changelog(prs: List):
    childprlist = []
    toremovelist = []
    for pr in prs:
        if has_label(pr, 'release notes: use body'):
            mdstring = body_to_release_notes(pr).strip()
            mdlines = mdstring.split('\n')
            pattern = r'\{.*\}$'
            for line in mdlines:
                cpr = copy.deepcopy(pr)
                mans = re.search(pattern, line)
                if mans:
                    label_list = mans.group().strip('{').strip('}').split(',')
                    for label in label_list:
                        label = label.strip()
                        if not (label in PRTYPES or label in TOPICS):
                            warning(
                                f"PR number #{pr['number']}'s changelog body has label {label}, "
                                "which is not a label we recognize ! We are ignoring this label. "
                                "This might result in a TODO changelog item!"
                            )
                            continue
                        cpr['labels'].append({'name': label})
                    mindex = mans.span()[0]
                    line = line[0:mindex]
                else:
                    warning(f"PR number #{pr['number']} is tagged as \"Use Body\", but the body "
                            "does not provide tags! This will result in TODO changelog items!")
                cpr['body'] = f'{line.strip()}\n'
                childprlist.append(cpr)
                if pr not in toremovelist:
                    toremovelist.append(pr)
    prs.extend(childprlist)
    prlist = [pr for pr in prs if pr not in toremovelist]
    return prlist

def main(new_version: str) -> None:
    major, minor, patchlevel = map(int, new_version.split("."))
    extra = ""
    release_type = 0 # 0 by default, 1 for point release, 2 for patch release
    if patchlevel == 0:
        # "major" release which changes just the minor version
        release_type = 1
        previous_minor = minor - 1
        basetag = f"v{major}.{minor}dev"
        # *exclude* PRs backported to previous stable-1.X branch
        #extra = f'-label:"backport {major}.{previous_minor}.x done"'
    else:
        # "minor" release which changes just the patchlevel
        release_type = 2
        previous_patchlevel = patchlevel - 1
        basetag = f"v{major}.{minor}.{previous_patchlevel}"
        # *include* PRs backported to current stable-4.X branch
        #extra = f'label:"backport {major}.{minor}.x done"'

    if release_type == 2:
        timestamp = get_tag_date(basetag)
    else:
        # Find the timestamp of the last shared commit
        shared_commit = subprocess.run([
            "git",
            "merge-base",
            basetag,
            "HEAD"
        ], shell=False, check=True, capture_output=True).stdout.decode().strip()
        timestamp = subprocess.run([
            "git",
            "show",
            "-s",
            "--format=\"%cI\"",
            shared_commit
        ], shell=False, check=True, capture_output=True).stdout.decode().strip().replace('"', '')
    print("Base tag is", basetag)
    print("Last common commit at ", timestamp)

    print("Downloading filtered PR list")
    prs = get_pr_list(timestamp, extra)
    prs = split_pr_into_changelog(prs)
    # print(json.dumps(prs, sort_keys=True, indent=4))

    # reset changelog file to state tracked in git

    subprocess.run(f'git checkout -- {FINALFILE}'.split(), check=True)

    changes_overview(prs, new_version)


if __name__ == "__main__":
    # the argument is the new version
    if len(sys.argv) == 1:
        itag = subprocess.run(
            [
                "gh",
                "release",
                "list",
                "--json=name,isLatest",
                "-q",
                ".[] | select(.isLatest == true)"
            ],
            shell=False,
            check=True,
            capture_output=True
        )
        itag = json.loads(itag.stdout.decode())["name"][1:]
        itag = itag.split('.')
        itag[-1] = str(int(itag[-1])+1)
        itag = ".".join(itag)
        main(itag)
    elif len(sys.argv) != 2:
        usage(sys.argv[0])
    else:
        main(sys.argv[1])
