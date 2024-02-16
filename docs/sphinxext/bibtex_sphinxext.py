"""Exposes a bibtex file to sphinx using bibtex-js under the covers.
This is largely based on the Jupyter notebook sphinx extension, that
you can find at https://github.com/jupyter/nbconvert

:Author: Anthony Scopatz <scopatz@gmail.com>
"""
import re
import os.path
import shutil
from codecs import open
from sphinx.util import ensuredir

from docutils import nodes
from docutils import utils
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives import unchanged_required, unchanged, flag


# Style is based off off SIAM
BASE_STYLE = """<div class="bibtex_template">
<li>
  <span class="if author">
    <span class="author"></span>.
  </span>
  "<span class="title"></span>,"
  <span style="if journal">
    <span class="journal" style="font-style: italic;"></span>
  </span>
  <span style="if booktitle">
    <span class="booktitle" style="font-style: italic;"></span>,
  </span>
  <span style="if school">
    <span class="school"></span>
  </span>
  <span class="if volume">
    <span class="volume"></span>
    <span class="if number">
      (<span class="number"></span>)
    </span>
    <span class="if pages">
        :<span class="pages"></span>
    </span>
  </span>
  <span class="if edition">
    <span class="edition"></span> ed.,
  </span>
  <span class="if month">
    <span class="month"></span>.
  </span>
  <span class="if year">
    <span class="year"></span>.
  </span>
  <span class="if url">
    <a class="url" style="color:black;">(view online)</a>
  </span>
</li>
</div>
"""

ENTRY_TYPES = [
    "article",
    "book",
    "booklet",
    "conference",
    "inbook",
    "incollection",
    "inproceedings",
    "manual",
    "mastersthesis",
    "misc",
    "phdthesis",
    "proceedings",
    "techreport",
    "unpublished",
]

nested_brace_pattern = re.compile("({.*?){(.*?)}(.*?})")
nested_brace_replace = lambda m: "".join(m.groups())


class Bibtex(Directive):
    """Uses bibtex-js to expose a bibtex file to sphinx.  This directive has
    one required argument, the path to the bibtex file.  It also may optionally
    have content that is raw html indicating the reference style.  For example,
    here is an "author-only" style::

        .. bibtex:: ../path/to/refs.bib

            <div class="bibtex_template">
            <div class="if author" style="font-weight: bold;">
              <span class="author"></span>
            </div>
            </div>

    These styles are defined by the bibtex-js project and more information may
    be found on their website [1].  The following explanation of style is taken
    from there:

        The default style can be specified explicitly by adding the following
        to your html code:

        .. code-block:: html

            <div class="bibtex_template">
            <div class="if author" style="font-weight: bold;">
              <span class="if year">
                <span class="year"></span>,
              </span>
              <span class="author"></span>
              <span class="if url" style="margin-left: 20px">
                <a class="url" style="color:black; font-size:10px">(view online)</a>
              </span>
            </div>
            <div style="margin-left: 10px; margin-bottom:5px;">
              <span class="title"></span>
            </div>
            </div>

        When class if is listed, the html element is only displayed if all fields
        listed as classes are present for an entry. Thus, the classes listed in
        <span class="if url"> cause the url related elements to only be displayed
        if an entry has the url field. For all other elements, the contents is
        replaced by the field-name specified as its class.

    This directive also has a 'no-count' options which may be used to supress
    the number of publication entries at the top::

        .. bibtex:: ../path/to/refs.bib
            :no-count:

    1. http://code.google.com/p/bibtex-js/wiki/styles
    """

    required_arguments = 1
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = {"no-count": flag}
    has_content = True

    def run(self):
        # check if raw html is supported
        if not self.state.document.settings.raw_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)

        # set up encoding
        attributes = {"format": "html"}
        show_count = "no-count" not in self.options

        # get path to bibtex file
        source_file = self.state.document.current_source
        source_dir = os.path.dirname(os.path.abspath(source_file))
        bibtex_path = os.path.normpath(os.path.join(source_dir, self.arguments[0]))
        bibtex_path = utils.relative_path(None, bibtex_path)

        # get and sanitize the bibtex
        with open(bibtex_path, encoding="utf-8") as f:
            bibtext = f.read()
        bibtext = bibtext.replace(r"\&", "&")
        bibtext = bibtext.replace(r"\bf", "")
        n = 1
        while 0 < n:
            bibtext, n = nested_brace_pattern.subn(nested_brace_replace, bibtext)
        if show_count:
            lower_bibtext = bibtext.lower()
            num_entries = 0
            for entry_type in ENTRY_TYPES:
                num_entries += lower_bibtext.count("@" + entry_type)

        styletext = BASE_STYLE if 0 == len(self.content) else "\n".join(self.content)

        text = '<textarea id="bibtex_input" style="display:none;">\n'
        text += bibtext
        text += "</textarea>\n"
        text += styletext
        if show_count:
            text += (
                '<a style="font-style: bold;">Number of entries: '
                "{0}</a><br/><br/>\n".format(num_entries)
            )
        text += '<div id="bibtex_display"></div>\n'

        # add dependency
        self.state.document.settings.record_dependencies.add(bibtex_path)
        attributes["source"] = bibtex_path

        # create notebook node
        bibtex_node = bibtex("", text, **attributes)
        (bibtex_node.source, bibtex_node.line) = self.state_machine.get_source_and_line(
            self.lineno
        )

        return [bibtex_node]


class bibtex(nodes.raw):
    pass


def visit_bibtex_node(self, node):
    self.visit_raw(node)


def depart_bibtex_node(self, node):
    self.depart_raw(node)


def setup(app):
    app.add_node(bibtex, html=(visit_bibtex_node, depart_bibtex_node))
    app.add_directive("bibtex", Bibtex)
    # Sphinx already comes with jquery, otherwise the following is needed.
    # app.add_javascript('http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/'
    #                   'jquery.min.js')
    bibjs = "bibtex_js.js"
    app.add_javascript(bibjs)
    source = os.path.abspath(os.path.join(os.path.dirname(__file__), bibjs))
    target = os.path.join(app.outdir, "_static", bibjs)
    ensuredir(os.path.dirname(target))
    shutil.copyfile(source, target)
