.. _theorymanual_bins:

===============================
Binning 
===============================

.. currentmodule:: pyne.tbins

This page explains the mathematics behind the binning functionality found in 
:py:mod:`bins`. The binning functions are purely mathematical in nature though 
they do have application in other parts of pyne that deal more directly physics 
with physics.

*****************************
Pointwise Linear Collapse
*****************************

The :func:`pointwise_linear_collapse` function takes an array of pointwise 
data :math:`y` that has the independent variable :math:`x` and collapses this into 
:math:`G` groups as defined by the bin boundaries :math:`x_g` where :math:`g` 
indexes :math:`G`. Both :math:`x` and :math:`x_g` must be monotonic and in the 
same direction.  Say that there are :math:`N` points in :math:`x` and :math:`y`.
Let :math:`n` index :math:`N`. Then for all points internal to a group, the 
collapsed value :math:`y_g` is: 

.. math::

    y_g = \frac{1}{x_{g+1} - x_g} \Sum_{n|x_g \le x_n}^{x_n \le x_{g+1}}
        \frac{y_{n+1} + y_n}{2} * (x_{n+1} - x_n)

The term :math:`(y_{n+1} + y_n)/2` is the center (average) value of a linear 
interpolation between the two points.  Therefore, :math:`y_g` is the 
:math:`x`-weighted average of :math:`y` over the entire group.  

In the event that the line between :math:`y_n` and :math:`y_{n+1}` crosses 
either the lower or upper bin boundary (or both) then their values are 
adjusted via a linear interpolation to the value at the bin boundary.

For a lower boundary crossing, the following substitutions are made to the 
equation above: 

.. math::

    x_n \to x_g

.. math::

    y_n \to \frac{y_{n+1} - y_n}{x_{n+1} - x_n} (x_g - x_n) + y_n

For an upper boundary crossing: 

.. math::

    x_{n+1} \to x_{g+1}

.. math::

    y_{n+1} \to \frac{y_{n+1} - y_n}{x_{n+1} - x_n} (x_{g+1} - x_n) + y_n

