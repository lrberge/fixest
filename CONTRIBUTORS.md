`fixest` is the outcome of many many years of development, but development alone in a cave. 
Github made it a collaborative project, channeling the attention of dozens of persons who took the time to report bugs, suggest new features and even make direct additions to the code.
The overall software quality is now nowhere near what it would have been without this global collaboration.

This document reports direct and indirect contributors to the project whom Laurent Berge (and likely the community) wishes to thank.<sup>[1](#fnt1)</sup>

## Direct contributors (PRs):

 - Sebastian Krantz: PR#66 substantial improvement of the `demean` function (hunted down overheads making it efficient for small data sizes, improved NA handling, fixed bugs and added user-friendly features).
 
 - Kevin Wilson: PR#78 automated the benchmarks. PR#79 implemented automatic checks via Github actions.
 
 - Grant McDermott: PR#130 rewrote substantially the introductory vignette.
 
## Indirect contributors:

Through issue filing, countless bugs have been resolved and many features have been added. In alphabetical order, a tentative and incomplete list of persons to thank: @adamaltmejd, @adamtheising, @amarbler, @apoorvalal, @benzipperer, @bgchamps, @chenwang, @clerousset, @clukewatson, @colejharvey, @d712, @dlindzee, @edrubin, @fostermeijer, @joseph-richard-martinez, @jurojas5, @kdzhang, @kendonB, @lyifa, @marissachilds, @nikolassch, @noahmbuckley, @nreigl, @Oravishayrizi, @pbaylis, @pei-huang, @phisherblack, @poliquin, @reifjulian, @rrichmond, @seunghoon001, @shoonlee, @SuperMayo, @tcovert, @tholdaway, @waynelapierre and @zozotintin. 

Special thanks to @karldw for writing the `broom` method very early, to @zozotintin for helping to debug problems occurring at the billion observations scale, and of course to Grant McDermott for single-handedly popularizing the package.

----

<a name="fnt1">1</a>: The distinction between direct and indirect contributors is quite arbitrary and does not pay tribute to some issue-filers who spent a lot of time and effort to resolve problems. But the rule, having made a PR, has the merit to be clear and simple.


