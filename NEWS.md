What’s new in version 1.2.1
---------------------------

### New functions added for root and extrema estimation

-   `classify_curve()` classifies a curve by convexity and shape type
-   `findmaxbell()` finds maximum for a bell curve
-   `findmaxtulip()` finds maximum for a tulip curve
-   `symextreme()` finds maximum for a symmetric curve
-   `findextreme()` finds extreme point for a general curve
-   `findroot()` finds root of a curve
-   `scan_curve()` scans a not noisy curve and finds all roots and
    extrema, inflections between them
-   `scan_noisy_curve()` scans a noisy curve and finds all roots and
    extrema, inflections between them

It is the implementation of Bell Extreme Finding Estimator (BEFE),
Integration Extreme Finding Estimator (IEFE) and Integration Root
Finding Estimator (IRFE) for roots, extrema and inflections of a curve,
see:

-   Demetris T. Christopoulos (2019). New methods for computing extremes
    and roots of a planar curve: introducing Noisy Numerical Analysis
    (2019). ResearchGate.
    <a href="http://dx.doi.org/10.13140/RG.2.2.17158.32324" class="uri">http://dx.doi.org/10.13140/RG.2.2.17158.32324</a>
