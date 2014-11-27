======================================
PyNE - The Nuclear Engineering Toolkit
======================================

.. raw:: html

    <div style="text-align:center;">

.. container:: frontpage-images

    .. image:: gallery/data_sources_thumb.png
    .. image:: gallery/fng_model_thumb.png
    .. image:: gallery/discretized_teapot_thumb.png
    .. image:: gallery/ace_thumb.png

.. raw:: html

    </div>

PyNE is a suite of tools to aid in
computational nuclear science & engineering.  PyNE seeks to provide
native implementations of common nuclear algorithms, as well as Python
bindings and I/O support for other industry standard nuclear codes.

.. raw:: html

    <script language="javascript"> 
    function pyneToggle(title, showHideDiv, switchTextDiv) {
        var ele = document.getElementById(showHideDiv);
        var text = document.getElementById(switchTextDiv);
        if(ele.style.display == "block") {
                ele.style.display = "none";
            text.innerHTML = title + " [+]";
        }
        else {
            ele.style.display = "block";
            text.innerHTML = title + " [-]";
        }
    } 
    </script>

    <div id="pynemenuheader">
      <a id="startHeader" 
         href="javascript:pyneToggle('Getting Started', 'startContent', 'startHeader');">
         Getting Started [+]</a>
    </div>
    <div style="clear:both;"></div>
    <div id="pynemenucontent">
      <div id="startContent" style="display:none;">
        <ul>
          <li><a href="install/index.html">Install</a></li>
          <li><a href="tutorial/index.html">Tutorial</a></li>
        </ul>
      </div>
    </div>

    <br />
    <div id="pynemenuheader">
      <a id="usingHeader" 
         href="javascript:pyneToggle('Using PyNE', 'usingContent', 'usingHeader');">
         Using PyNE [+]</a>
    </div>
    <div style="clear:both;"></div>
    <div id="pynemenucontent">
      <div id="usingContent" style="display:none;">
        <ul>
          <li><a href="usersguide/index.html">User's Guide</a></li>
          <li><a href="pyapi/index.html}}">Python API Documentation</a></li>
          <li><a href="cppapi/index.html">C++ & Fortran API Documentation</a></li>
          <li><a href="mailto:pyne-users+subscribe@googlegroups.com?subject=Subscribe&body=Send this message to subscribe to the list">Join</a> the <a href="https://groups.google.com/forum/#!forum/pyne-users" target="_blank"> Users</a> mailing list.
          <li><a href="https://github.com/pyne/pyne/issues">Report an Issue</a></li>
        </ul>
      </div>
    </div>

    <br />
    <div id="pynemenuheader">
      <a id="contributeHeader" 
         href="javascript:pyneToggle('Contribute', 'contributeContent', 'contributeHeader');">
         Contribute [+]</a>
    </div>
    <div style="clear:both;"></div>
    <div id="pynemenucontent">
      <div id="contributeContent" style="display:none;">
        <ul>
          <li><a href="devsguide/index.html">Developer's Guide</a></li>
          <li><a href="http://github.com/pyne/pyne">Source Code</a></li>
          <li><a href="mailto:pyne-users+subscribe@googlegroups.com?subject=Subscribe&body=Send this message to subscribe to the list">Join</a> the 
              <a href="https://groups.google.com/forum/#!forum/pyne-users" target="_blank">Developers</a> mailing list.
          <li><a href="dev_team.html">The PyNE Team</a></li>
        </ul>
      </div>
    </div>

    <br />
    <div id="pynemenuheader">
      <a id="learnHeader" 
         href="javascript:pyneToggle('Learn More', 'learnContent', 'learnHeader');">
         Learn More [+]</a>
    </div>
    <div style="clear:both;"></div>
    <div id="pynemenucontent">
      <div id="learnContent" style="display:none;">
        <ul>
          <li><a href="theorymanual/index.html">Theory Manual</a></li>
          <li><a href="pubs.html">Publications</a></li>
          <li><a href="previous/index.html">Release Notes</a></li>
          <li><a href="gsoc/index.html">Project Ideas</a></li>
        </ul>
      </div>
    </div>

.. toctree::
     :maxdepth: 1
     install
     tutorial
     users_guide
     devsguide
     theorymanual
     pyapi
     cppapi
..  
..      gallery/index
..      previous/index
..      dev_team
..      pubs
..      gsoc/index
..  
.. _C++ API: cppapi/html/index.html

.. _GitHub project site: https://github.com/pyne

.. _github: https://github.com/pyne/pyne

.. _zip: https://github.com/pyne/pyne/zipball/0.4
.. _tar: https://github.com/pyne/pyne/tarball/0.4
