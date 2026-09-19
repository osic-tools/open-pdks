#!/usr/bin/env python3
#
# create_lib_library.py
#
#----------------------------------------------------------------------------
# Given a destination directory holding individual liberty files of a number
# of cells, create a single liberty library file named <alllibname> and place
# it in the same directory.  This is done for the option "compile" if specified
# for the "-lib" install.
#
# Because liberty files tend to be both huge and ASCII, this script supports
# handling sets of gzipped .lib files (expecting name ".lib.gz").
#
# The individual files may be either complete liberty files (each with its
# own "library (...) { ... }" header wrapped around the cell definitions),
# or contain only cell definitions.  If the files contain only cell
# definitions, the header must exist in its own file and be passed to this
# script with the "-header" option.  When files do have their own headers,
# they are merged together so that, e.g., any templates that are defined in
# some files but not in others will appear in the combined library header.
#
# Original python script by Tim Edwards for the original open_pdks.
#
# Updated September 17, 2026
# Modifications/corrections by Tim Edwards and additional parsing routines
# by Claude Opus to implement a handful of wish-list items such as merging
# multiple headers together into a single library header.  The new code also
# corrects an issue that the original script had with failing to remove a
# closing brace from the end of individual self-contained liberty files.
#----------------------------------------------------------------------------

import sys
import os
import re
import gzip
import glob
import fnmatch
import natural_sort

#----------------------------------------------------------------------------

def usage():
    print('')
    print('Usage:')
    print('    create_lib_library <destlibdir> <destlib> [-compile-only] ')
    print('             [-excludelist="file1,file2,..."]')
    print('')
    print('Create a single liberty library from a set of individual liberty files.')
    print('')
    print('where:')
    print('    <destlibdir>      is the directory containing the individual liberty files')
    print('    <destlib>         is the root name of the library file')
    print('    -compile-only     remove the indidual files if specified')
    print('    -excludelist=     is a comma-separated list of files to ignore')
    print('    -header=		 is the name of a file containing header information')
    print('')

#----------------------------------------------------------------------------
# Warning:  This script is unfinished.  The library name in the header needs
# to be changed to the full library name.  Also:  There is no mechanism for
# collecting all files belonging to a single process corner/temperature/
# voltage.
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Minimal liberty tokenizer.  This does not attempt to understand liberty
# syntax, which comes in many valid styles.  It only recognizes quoted
# strings, comments, and line continuations (so that punctuation inside them
# is ignored) and the punctuation needed to find statement and group
# boundaries.  Everything else is passed through as "word" tokens.
#----------------------------------------------------------------------------

token_re = re.compile(r'''
      (?P<string>"(?:\\.|[^"\\])*")
    | (?P<comment>/\*.*?\*/|//[^\r\n]*)
    | (?P<cont>\\[ \t]*\r?\n)
    | (?P<newline>\r?\n)
    | (?P<space>[ \t\r\f\v]+)
    | (?P<punct>[{}();])
    | (?P<word>[^\s"{}();\\/]+|[\\/])
''', re.S | re.X)

brace_re = re.compile(r'"(?:\\.|[^"\\])*"|/\*.*?\*/|//[^\r\n]*|[{}]', re.S)

# C-style "//" comments are not part of the liberty format.  Some tools
# (e.g., yosys) accept them, but others (e.g., OpenSTA) do not, so they are
# converted to C-style comments in the output.

line_comment_re = re.compile(r'"(?:\\.|[^"\\])*"|/\*.*?\*/|//([^\r\n]*)', re.S)

number_re = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?')

ignored_tokens = ('comment', 'space', 'cont')

# Library-level complex attributes that may legitimately appear more than
# once, and how to tell the instances apart.  All other attributes are
# expected to appear only once, and a different value in another file's
# header is reported as a conflict.

multi_firstarg = ('voltage_map',)
multi_allargs = ('define', 'define_group', 'define_cell_area')

#----------------------------------------------------------------------------
# Parse the header of a liberty file.
#
# Returns None if the text does not begin with a "library" group (such as a
# file containing only cell definitions).  Otherwise returns a dictionary:
#
#   'cell_start':  Index of the first cell group, or of the library group's
#                  closing brace if there are no cells.
#   'complete':    True if the first cell or the end of the library group
#                  was found (False if the text ended first).
#   'at_close':    True if 'cell_start' is the library group's closing brace.
#   'statements':  List of (name, args, is_group, start, end) for each
#                  library-level statement preceding the first cell.
#----------------------------------------------------------------------------

def parse_liberty_header(text):
    tokens = token_re.finditer(text)

    # The first significant token must be the "library" keyword.
    for tok in tokens:
        if tok.lastgroup in ignored_tokens or tok.lastgroup == 'newline':
            continue
        if tok.lastgroup == 'word' and tok.group() == 'library':
            break
        return None
    else:
        return None

    # Find the opening brace of the library group
    for tok in tokens:
        if tok.lastgroup == 'punct':
            if tok.group() == '{':
                break
            elif tok.group() in ';}':
                return None
    else:
        return None

    result = {'cell_start': len(text), 'complete': False, 'at_close': False,
		'statements': []}
    statements = result['statements']

    start = None	# Start of the current statement
    name = None		# Name of the current statement
    args = None		# Text inside the statement's first parentheses
    argstart = None
    paren = 0
    pending = None	# Possible end of a statement without a semicolon

    for tok in tokens:
        kind = tok.lastgroup
        value = tok.group()
        if kind in ignored_tokens:
            continue
        if kind == 'newline':
            if start is not None and paren == 0 and pending is None:
                pending = tok.start()
            continue

        # A statement not terminated by a semicolon ends at the newline,
        # unless the next line begins with the opening brace of a group.
        if pending is not None:
            if kind == 'punct' and value == '{':
                pending = None
            else:
                statements.append((name, args, False, start, pending))
                start = None
                pending = None

        if start is None:
            if kind == 'punct':
                if value == '}':
                    result['cell_start'] = tok.start()
                    result['complete'] = True
                    result['at_close'] = True
                    return result
                # Ignore stray punctuation between statements
                continue
            if value == 'cell' or value == 'scaled_cell':
                result['cell_start'] = tok.start()
                result['complete'] = True
                return result
            start = tok.start()
            name = value
            args = None
            argstart = None
            paren = 0
            continue

        if kind != 'punct':
            continue

        if value == '(':
            paren += 1
            if paren == 1 and argstart is None:
                argstart = tok.end()
        elif value == ')':
            if paren > 0:
                paren -= 1
                if paren == 0 and args is None:
                    args = text[argstart:tok.start()]
        elif value == ';' and paren == 0:
            statements.append((name, args, False, start, tok.end()))
            start = None
        elif value == '{' and paren == 0:
            depth = 1
            for tok in tokens:
                if tok.lastgroup == 'punct':
                    if tok.group() == '{':
                        depth += 1
                    elif tok.group() == '}':
                        depth -= 1
                        if depth == 0:
                            break
            else:
                # Group was not closed before the end of the text
                return result
            statements.append((name, args, True, start, tok.end()))
            start = None
        elif value == '}':
            # Closing brace of the library group ends an unterminated statement
            statements.append((name, args, False, start, tok.start()))
            result['cell_start'] = tok.start()
            result['complete'] = True
            result['at_close'] = True
            return result

    return result

#----------------------------------------------------------------------------
# Find the closing brace of the library group, given the index of a
# position inside the library group at brace depth 1 (such as the start of
# the first cell).  Returns None if the group is not closed.
#----------------------------------------------------------------------------

def find_library_close(text, pos):
    depth = 1
    for match in brace_re.finditer(text, pos):
        value = match.group()
        if value == '{':
            depth += 1
        elif value == '}':
            depth -= 1
            if depth == 0:
                return match.start()
    return None

#----------------------------------------------------------------------------
# Normalize a statement for comparison, so that differences in whitespace,
# comments, line continuations, optional semicolons, quoting, and number
# formatting are ignored.
#----------------------------------------------------------------------------

def normalize_statement(text):
    words = []
    for tok in token_re.finditer(text):
        if tok.lastgroup in ignored_tokens or tok.lastgroup == 'newline':
            continue
        value = tok.group()
        if value == ';' or value == '\\':
            continue
        if tok.lastgroup == 'string':
            value = value[1:-1]
        words.append(value)
    ntext = ' '.join(words)
    ntext = number_re.sub(lambda m: repr(float(m.group())), ntext)
    ntext = re.sub(r'\s*([,:;(){}])\s*', r'\1', ntext)
    ntext = re.sub(r'\s+', ' ', ntext).strip()
    return ntext.rstrip(';')

#----------------------------------------------------------------------------
# Key used to decide whether two library-level statements from different
# files define the same thing.
#----------------------------------------------------------------------------

def statement_key(name, args, is_group):
    if is_group:
        return (name, normalize_statement(args or ''))
    elif args is not None and name in multi_allargs:
        return (name, normalize_statement(args))
    elif args is not None and name in multi_firstarg:
        return (name, normalize_statement(args.split(',')[0]))
    else:
        return (name,)

#----------------------------------------------------------------------------
# Return the text of a statement, including any indentation preceding it on
# the same line.
#----------------------------------------------------------------------------

def statement_text(text, start, end):
    linestart = text.rfind('\n', 0, start) + 1
    if text[linestart:start].strip() == '':
        start = linestart
    return text[start:end]

#----------------------------------------------------------------------------
# Convert any "//" comments in liberty text to "/* */" comments.
#
# This probably should raise a warning or error rather than be quietly
# corrected;  line-based comments are not part of the official liberty
# syntax.
#----------------------------------------------------------------------------

def convert_line_comments(text):
    if '//' not in text:
        return text

    def replace(match):
        if match.group(1) is None:
            # Quoted string or C-style comment;  leave unchanged
            return match.group(0)
        return '/*' + match.group(1).replace('*/', '* /') + ' */'

    return line_comment_re.sub(replace, text)

#----------------------------------------------------------------------------
# Read a liberty file.  If header_only is True, read only as much of the file
# as needed to parse the header (liberty files can be very large).
#----------------------------------------------------------------------------

def read_liberty(lfile, compressed, header_only=False):
    if compressed:
        ifile = gzip.open(lfile, 'rt')
    else:
        ifile = open(lfile, 'r')

    with ifile:
        if not header_only:
            return ifile.read()

        cell_re = re.compile(r'\b(scaled_)?cell\s*\(')
        chunks = []
        for line in ifile:
            chunks.append(line)
            if cell_re.search(line):
                # This may be a cell mentioned in a comment, in which case
                # parsing will be incomplete and reading continues.
                text = ''.join(chunks)
                parsed = parse_liberty_header(text)
                if parsed is None or parsed['complete']:
                    return text
        return ''.join(chunks)

#----------------------------------------------------------------------------
# Compile a liberty library file from constituent parts in separate files
#----------------------------------------------------------------------------

def create_lib_library(destlibdir, destlib, do_compile_only=False, excludelist=[],
	headerfile=None):

    compressed = False

    # destlib should not have a file extension
    destlibroot = os.path.splitext(destlib)[0]

    # If a header file is specified, read it first.  If the header file
    # name is the same as the library name (typical), the file will be
    # destroyed before being rebuilt, so save the header file contents.

    htext = None
    if headerfile:
        try:
            with open(destlibdir + '/' + headerfile, 'r') as ifile:
                htext = ifile.read()
        except:
            print('Error reading liberty header file ' + headerfile)
            headerfile = None

    alllibname = destlibdir + '/' + destlibroot + '.lib'
    if os.path.isfile(alllibname):
        os.remove(alllibname)

    print('Diagnostic:  Creating consolidated liberty library ' + destlibroot + '.lib')

    # If file "filelist.txt" exists in the directory, get the list of files from it
    if os.path.exists(destlibdir + '/filelist.txt'):
        with open(destlibdir + '/filelist.txt', 'r') as ifile:
            rlist = ifile.read().splitlines()
            llist = []
            for rfile in rlist:
                llist.append(destlibdir + '/' + rfile)
    else:
        llist = glob.glob(destlibdir + '/*.lib.gz')
        if len(llist) == 0:
            llist = glob.glob(destlibdir + '/*.lib')
        else:
            compressed = True
        llist = natural_sort.natural_sort(llist)

    # Create exclude list with glob-style matching using fnmatch
    if len(llist) > 0:
        llistnames = list(os.path.split(item)[1] for item in llist)
        notllist = []
        for exclude in excludelist:
            notllist.extend(fnmatch.filter(llistnames, exclude))

        # Apply exclude list
        if len(notllist) > 0:
            for file in llist[:]:
                if os.path.split(file)[1] in notllist:
                    llist.remove(file)

    if len(llist) > 1:
        print('New file is:  ' + alllibname)

        # Pass 1:  Build the library header.  The base header is taken from
        # the header file if one was given, or else from the first file that
        # has a header.  Library-level definitions (such as table templates)
        # found in the headers of the other files and not already present
        # are added to it.

        header = None		# Text of the base header
        headersource = None	# Where the base header came from
        headerkeys = {}		# Statement key -> (normalized text, source)
        extras = []		# Statements added from other files' headers
        conflicts = set()

        if htext is not None:
            parsed = parse_liberty_header(htext)
            if parsed is None:
                # Not parseable as a library group;  use it as-is, minus
                # the closing brace of the library group.
                print('Warning:  Could not parse header file ' + headerfile +
			';  header will not be merged with other files.')
                header = '\n'.join(line for line in htext.splitlines()
			if not line.startswith('}'))
            else:
                header = htext[:parsed['cell_start']]
                for (name, args, is_group, start, end) in parsed['statements']:
                    key = statement_key(name, args, is_group)
                    headerkeys[key] = (normalize_statement(htext[start:end]), headerfile)
            headersource = headerfile

        for lfile in llist:
            if not os.path.exists(lfile):
                continue
            ltext = read_liberty(lfile, compressed, header_only=True)
            parsed = parse_liberty_header(ltext)
            if parsed is None:
                continue
            lname = os.path.split(lfile)[1]

            if header is None:
                header = ltext[:parsed['cell_start']]
                headersource = lname
                for (name, args, is_group, start, end) in parsed['statements']:
                    key = statement_key(name, args, is_group)
                    headerkeys[key] = (normalize_statement(ltext[start:end]), lname)
                continue

            for (name, args, is_group, start, end) in parsed['statements']:
                key = statement_key(name, args, is_group)
                ntext = normalize_statement(ltext[start:end])
                if key not in headerkeys:
                    headerkeys[key] = (ntext, lname)
                    extras.append(statement_text(ltext, start, end))
                elif headerkeys[key][0] != ntext and key not in conflicts:
                    conflicts.add(key)
                    what = name if len(key) == 1 else name + ' (' + key[1] + ')'
                    print('Warning:  Header definition of ' + what + ' in ' + lname +
			' differs from the one in ' + headerkeys[key][1] +
			';  using the one in ' + headerkeys[key][1] + '.')

        if len(extras) > 0:
            print('Merged ' + str(len(extras)) + ' header definitions from other files into header from ' + headersource + '.')

        # Pass 2:  Write the header, then the cell definitions from each file.

        with open(alllibname, 'w') as ofile:

            if header is not None:
                print(convert_line_comments(header.rstrip()), file=ofile)
                if len(extras) > 0:
                    print('\n/* Header definitions merged from other files */', file=ofile)
                    for extra in extras:
                        print(convert_line_comments(extra), file=ofile)
                print('', file=ofile)

            for lfile in llist:
                if not os.path.exists(lfile):
                    print('Error: File ' + lfile + ' not found (skipping).')
                    continue
                ltext = read_liberty(lfile, compressed)

                parsed = parse_liberty_header(ltext)
                if parsed is None:
                    # File contains only cell definitions;  copy verbatim
                    print(convert_line_comments(ltext.rstrip()), file=ofile)
                else:
                    # Strip the header and the library group's closing brace
                    cellstart = parsed['cell_start']
                    if parsed['at_close']:
                        cellend = cellstart
                    else:
                        cellend = find_library_close(ltext, cellstart)
                        if cellend is None:
                            print('Warning:  Library group in ' + lfile + ' is not closed.')
                            cellend = len(ltext)
                    celltext = statement_text(ltext, cellstart, cellend)
                    if celltext.strip() != '':
                        print(convert_line_comments(celltext.rstrip()), file=ofile)
                print('/*--------EOF---------*/\n', file=ofile)

            # Close the library group
            if header is not None:
                print('}', file=ofile)
            else:
                print('Warning:  No library header found in files or given by -header.')

        if do_compile_only == True:
            print('Compile-only:  Removing individual liberty files')
            for lfile in llist:
                if os.path.isfile(lfile):
                    os.remove(lfile)
    else:
        print('Only one file (' + str(llist) + ');  ignoring "compile" option.')

#----------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) == 1:
        usage()
        sys.exit(0)

    argumentlist = []

    # Defaults
    do_compile_only = False
    excludelist = []
    headerfile = None

    # Break arguments into groups where the first word begins with "-".
    # All following words not beginning with "-" are appended to the
    # same list (optionlist).  Then each optionlist is processed.
    # Note that the first entry in optionlist has the '-' removed.

    for option in sys.argv[1:]:
        if option.find('-', 0) == 0:
            keyval = option[1:].split('=')
            if keyval[0] == 'compile-only':
                if len(keyval) > 1:
                    if keyval[1].lower() == 'true' or keyval[1].lower() == 'yes' or keyval[1] == '1':
                        do_compile_only = True
                else:
                    do_compile_only = True
            elif keyval[0] == 'exclude' or keyval[0] == 'excludelist':
                if len(keyval) > 1:
                    excludelist = keyval[1].strip('"').split(',')
                else:
                    print("No items in exclude list (ignoring).")
            elif keyval[0] == 'header':
                if len(keyval) > 1:
                    headerfile = keyval[1].strip('"')
                else:
                    print("No value for header file (ignoring).")
            else:
                print("Unknown option '" + keyval[0] + "' (ignoring).")
        else:
            argumentlist.append(option)

    if len(argumentlist) < 2:
        print("Not enough arguments given to create_lib_library.py.")
        usage()
        sys.exit(1)

    destlibdir = argumentlist[0]
    destlib = argumentlist[1]

    print('')
    print('Create liberty library from files:')
    print('')
    print('Path to files: ' + destlibdir)
    print('Name of compiled library: ' + destlib + '.lib')
    print('Remove individual files: ' + ('Yes' if do_compile_only else 'No'))
    if len(excludelist) > 0:
        print('List of files to exclude: ')
        for file in excludelist:
            print(file)
    print('')

    create_lib_library(destlibdir, destlib, do_compile_only, excludelist, headerfile)
    print('Done.')
    sys.exit(0)

#----------------------------------------------------------------------------
