


print.dsb = function(x, ...){
    cat(x, sep = "\n")
}



#' Simple and powerful string manipulation with the dot square bracket operator
#'
#' Compactly performs many low level string operations. Advanced support for pluralization.
#'
#'
#'
#' @param ... Character scalars that will be collapsed with the argument \code{sep}. You can use \code{".[x]"} within each character string to insert the value of \code{x} in the string. You can add string operations in each \code{".[]"} instance with the syntax \code{"'arg'op ? x"} (resp. \code{"'arg'op ! x"}) to apply the operation \code{'op'} with the argument \code{'arg'} to \code{x} (resp. the verbatim of \code{x}). Otherwise, what to say? Ah, nesting is enabled, and since there's over 30 operators, it's a bit complicated to sort you out in this small space. But type \code{dsb("--help")} to prompt an (almost) extensive help.
#' @param frame An environment used to evaluate the variables in \code{".[]"}.
#' @param sep Character scalar, default is \code{""}. It is used to collapse all the elements in \code{...}.
#' @param vectorize Logical, default is \code{FALSE}. If \code{TRUE}, Further, elements in \code{...} are NOT collapsed together, but instead vectorised.
#' @param nest Logical, default is \code{TRUE}. Whether the original character strings should be nested into a \code{".[]"}. If \code{TRUE}, then things like \code{dsb("S!one, two")} are equivalent to \code{dsb(".[S!one, two]")} and hence create the vector \code{c("one", "two")}.
#' @param collapse Character scalar or \code{NULL} (default). If provided, the resulting character vector will be collapsed into a character scalar using this value as a separator.
#'
#'
#' There are over 30 basic string operations, it supports pluralization, it's fast (e.g. faster than \code{glue} in the benchmarks), string operations can be nested (it may be the most powerful feature), operators have sensible defaults.
#'
#' See detailed help on the console with \code{dsb("--help")}. The real help is in fact in the "Examples" section.
#'
#'
#' @return
#' It returns a character vector whose length depends on the elements and operations in \code{".[]"}.
#'
#' @examples
#'
#' #
#' # BASIC USAGE ####
#' #
#'
#' x = c("Romeo", "Juliet")
#'
#' # .[x] inserts x
#' dsb("Hello .[x]!")
#'
#' # elements in ... are collapsed with "" (default)
#' dsb("Hello .[x[1]], ",
#'     "how is .[x[2]] doing?")
#'
#' # Splitting a comma separated string
#' # The mechanism is explained later
#' dsb("/J. Mills, David, Agnes, Dr Strong")
#'
#' # Nota: this is equivalent to (explained later)
#' dsb("', *'S !J. Mills, David, Agnes, Dr Strong")
#'
#' #
#' # Applying low level operations to strings
#' #
#'
#' # Two main syntax:
#'
#' # A) expression evaluation
#' # .[operation ? x]
#' #             | |
#' #             |  \-> the expression to be evaluated
#' #              \-> ? means that the expression will be evaluated
#'
#' # B) verbatim
#' # .[operation ! x]
#' #             | |
#' #             |  \-> the expression taken as verbatim (here ' x')
#' #              \-> ! means that the expression is taken as verbatim
#'
#' # operation: usually 'arg'op with op an operation code.
#'
#' # Example: splitting
#' x = "hello dear"
#' dsb(".[' 's ? x]")
#' # x is split by ' '
#'
#' dsb(".[' 's !hello dear]")
#' # 'hello dear' is split by ' '
#' # had we used ?, there would have been an error
#'
#' # By default, the string is nested in .[], so in that case no need to use .[]:
#' dsb("' 's ? x")
#' dsb("' 's !hello dear")
#'
#' # There are 35 string operators
#' # Operators usually have a default value
#' # Operations can be chained by separating them with a comma
#'
#' # Example: default of 's' is ' ' + chaining with collapse
#' dsb("s, ' my 'c!hello dear")
#'
#' #
#' # Nesting
#' #
#'
#' # .[operations ! s1.[expr]s2]
#' #              |    |
#' #              |     \-> expr will be evaluated then added to the string
#' #               \-> nesting requires verbatim evaluation: '!'
#'
#' dsb("The variables are: .[C!x.[1:4]].")
#'
#' # This one is a bit ugly but it shows triple nesting
#' dsb("The variables are: .[w, C!.[2* ! x.[1:4]].[S, 4** ! , _sq]].")
#'
#' #
#' # Splitting
#' #
#'
#' # s: split with fixed pattern, default is ' '
#' dsb("s !a b c")
#' dsb("' b 's !a b c")
#'
#' # S: split with regex pattern, default is ', *'
#' dsb("S !a, b, c")
#' dsb("'[[:punct:] ]'S !a! b; c")
#'
#' #
#' # Collapsing
#' #
#'
#' # c and C do the same, their default is different
#' # syntax: 's1||s2' with
#' # - s1 the string used for collapsing
#' # - s2 (optional) the string used for the last collapse
#'
#' # c: default is ' '
#' dsb("c?1:3")
#'
#' # C: default is ', || and '
#' dsb("C?1:3")
#'
#' dsb("', || or 'c?1:4")
#'
#' #
#' # Extraction
#' #
#'
#' # x: extracts the first pattern
#' # X: extracts all patterns
#' # syntax: 'pattern'x
#' # Default is '[[:alnum:]]+'
#'
#' x = "This years is... 2020"
#' dsb("x ? x")
#' dsb("X ? x")
#'
#' dsb("'\\d+'x ? x")
#'
#' #
#' # STRING FORMATTING ####
#' #
#'
#' #
#' # u, U: uppercase first/all letters
#'
#' # first letter
#' dsb("u!julia mills")
#'
#' # title case: split -> upper first letter -> collapse
#' dsb("s, u, c!julia mills")
#'
#' # upper all letters
#' dsb("U!julia mills")
#'
#' #
#' # L: lowercase
#'
#' dsb("L!JULIA MILLS")
#'
#' #
#' # q, Q: single or double quote
#'
#' dsb("S, q, C!Julia, David, Wilkins")
#' dsb("S, Q, C!Julia, David, Wilkins")
#'
#' #
#' # f, F: formats the string to fit the same length
#'
#'
#' score = c(-10, 2050)
#' nm = c("Wilkins", "David")
#' dsb("Monopoly scores:\n.['\n'c ! - .[f ? nm]: .[F ? score] US$]")
#'
#' # OK that example may have been a bit too complex,
#' # let's make it simple:
#'
#' dsb("Scores: .[f ? score]")
#' dsb("Names: .[F ? nm]")
#'
#' #
#' # w, W: reformat the white spaces
#' # w: suppresses trimming white spaces + normalizes successive white spaces
#' # W: same but also includes punctuation
#'
#' dsb("w ! The   white  spaces are now clean.  ")
#'
#' dsb("W ! I, really -- truly; love punctuation!!!")
#'
#' #
#' # %: applies sprintf formatting
#'
#' dsb("pi = .['.2f'% ? pi]")
#'
#' #
#' # a: appends text on each item
#' # syntax: 's1|s2'a, adds s1 at the beginning and s2 at the end of the string
#' # It accepts the special values :1:, :i:, :I:, :a:, :A:
#' # These values create enumerations (only one such value is accepted)
#'
#' # appending square brackets
#' dsb("'[|]'a, ' + 'c!x.[1:4]")
#'
#' # Enumerations
#' acad = dsb("/you like admin, you enjoy working on weekends, you really love emails")
#' dsb("Main reasons to pursue an academic career:\n .[':i:) 'a, C ? acad].")
#'
#' #
#' # A: same as 'a' but adds at the begging/end of the full string (not on the elements)
#' # special values: :n:, :N:, give the number of elements
#'
#' characters = dsb("/David, Wilkins, Dora, Agnes")
#' dsb("There are .[':N: characters: 'A, C ? characters].")
#'
#'
#' #
#' # stop: removes basic English stopwords
#' # the list is from the Snowball project: http://snowball.tartarus.org/algorithms/english/stop.txt
#'
#' dsb("stop, w!It is a tale told by an idiot, full of sound and fury, signifying nothing.")
#'
#' #
#' # k: keeps the first n characters
#' # syntax: nk: keeps the first n characters
#' #         'n|s'k: same + adds 's' at the end of shortened strings
#' #         'n||s'k: same but 's' counts in the n characters kept
#'
#' words = dsb("/short, constitutional")
#' dsb("5k ? words")
#'
#' dsb("'5|..'k ? words")
#'
#' dsb("'5||..'k ? words")
#'
#' #
#' # K: keeps the first n elements
#' # syntax: nK: keeps the first n elements
#' #         'n|s'K: same + adds the element 's' at the end
#' #         'n||s'K: same but 's' counts in the n elements kept
#' #
#' # Special values :rest: and :REST:, give the number of items dropped
#'
#' bx = dsb("/Pessac Leognan, Saint Emilion, Marguaux, Saint Julien, Pauillac")
#' dsb("Bordeaux wines I like: .[3K, ', 'C ? bx].")
#'
#' dsb("Bordeaux wines I like: .['3|etc..'K, ', 'C ? bx].")
#'
#' dsb("Bordeaux wines I like: .['3||etc..'K, ', 'C ? bx].")
#'
#' dsb("Bordeaux wines I like: .['3|and at least :REST: others'K, ', 'C ? bx].")
#'
#' #
#' # Ko, KO: special operator which keeps the first n elements and adds "others"
#' # syntax: nKo
#' # KO gives the rest in letters
#'
#' dsb("Bordeaux wines I like: .[4KO, C ? bx].")
#'
#' #
#' # r, R: string replacement
#' # syntax: 's'R: deletes the content in 's' (replaces with the empty string)
#' #         's1 => s2'R replaces s1 into s2
#' # r: fixed / R: perl = TRUE
#'
#' dsb("'e'r !The letter e is deleted")
#'
#' # adding a perl look-behind
#' dsb("'(?<! )e'R !The letter e is deleted")
#'
#' dsb("'e => a'r !The letter e becomes a")
#'
#' dsb("'([[:alpha:]]{3})[[:alpha:]]+ => \\1.'R !Trimming the words")
#'
#' #
#' # *, *c, **, **c: replication, replication + collapse
#' # syntax: n* or n*c
#' # ** is the same as * but uses "each" in the replication
#'
#' dsb("N.[10*c!o]!")
#'
#' dsb("3*c ? 1:3")
#' dsb("3**c ? 1:3")
#'
#' #
#' # d: replaces the items by the empty string
#' # -> useful in conditions
#'
#' dsb("d!I am going to be annihilated")
#'
#' #
#' # ELEMENT MANIPULATION ####
#' #
#'
#' #
#' # D: deletes all elements
#' # -> useful in conditions
#'
#' x = dsb("/I'll, be, deleted")
#' dsb("D ? x")
#'
#' #
#' # i, I: inserts an item
#' # syntax: 's1|s2'i: inserts s1 first and s2 last
#' # I: is the same as i but is 'invisibly' included
#'
#' characters = dsb("/David, Wilkins, Dora, Agnes, Trotwood")
#' dsb("'Heep|Spenlow'i, C ? characters")
#'
#' dsb("'Heep|Spenlow'I, C ? characters")
#'
#'
#' #
#' # PLURALIZATION ####
#' #
#'
#' # There is support for pluralization
#'
#' #
#' # *s, *s_: adds 's' or 's ' depending on the number of elements
#'
#' nb = 1:5
#' dsb("Number.[*s, D ? nb]: .[C ? nb]")
#' dsb("Number.[*s, D ? 2 ]: .[C ? 2 ]")
#'
#' # or
#' dsb("Number.[*s, ': 'A, C ? nb]")
#'
#'
#' #
#' # v, V: adds a verb at the beginning/end of the string
#' # syntax: 'verb'v
#'
#' # Unpopular opinion?
#' brand = c("Apple", "Samsung")
#' dsb(".[V, C ? brand] overrated.")
#' dsb(".[V, C ? brand[1]] overrated.")
#'
#' win = dsb("/Peggoty, Agnes, Emily")
#' dsb("The winner.[*s_, v, C ? win].")
#' dsb("The winner.[*s_, v, C ? win[1]].")
#'
#' # Other verbs
#' dsb(".[' have'V, C ? win] won a prize.")
#' dsb(".[' have'V, C ? win[1]] won a prize.")
#'
#' dsb(".[' was'V, C ? win] unable to come.")
#' dsb(".[' was'V, C ? win[1]] unable to come.")
#'
#' #
#' # *A: appends text depending on the length of the vector
#' # syntax: 's1|s2 / s3|s4'
#' #         if length == 1: applies 's1|s2'A
#' #         if length >  1: applies 's3|s4'A
#'
#' win = dsb("/Barkis, Micawber, Murdstone")
#' dsb("The winner.[' is /s are '*A, C ? win].")
#' dsb("The winner.[' is /s are '*A, C ? win[1]].")
#'
#' #
#' # CONDITIONS ####
#' #
#'
#' # Conditions can be applied with 'if' statements.",
#' # The syntax is 'type comp value'if(true : false), with
#' # - type: either 'len', 'char', 'fixed' or 'regex'
#' #   + len: number of elements in the vector
#' #   + char: number of characters
#' #   + fixed: fixed pattern
#' #   + regex: regular expression pattern
#' # - comp: a comparator:
#' #   + valid for len/char: >, <, >=, <=, !=, ==
#' #   + valid for fixed/regex: !=, ==
#' # - value: a value for which the comparison is applied.
#' # - true: operations to be applied if true (can be void)
#' # - false: operations to be applied if false (can be void)
#'
#' dsb("'char <= 2'if('(|)'a : '[|]'a), ' + 'c ? c(1, 12, 123)")
#'
#' sentence = "This is a sentence with some longish words."
#' dsb("s, 'char<=4'if(D), c ? sentence")
#'
#' dsb("s, 'fixed == e'if(:D), c ! Only words with an e are selected.")
#'
#' #
#' # ARGUMENTS FROM THE FRAME ####
#' #
#'
#' # Arguments can be evaluated from the calling frame.
#' # Simply use backticks instead of quotes.
#'
#' dollar = 6
#' reason = "glory"
#' dsb("Why do you develop packages? For .[`dollar`*c!$]?",
#'     "For money? No... for .[U,''s, c?reason]!", sep = "\n")
#'
#'
#'
#'
#'
dsb = function(..., frame = parent.frame(), sep = "", vectorize = FALSE, nest = TRUE,
               collapse = NULL){
    check_arg(vectorize, nest, "logical scalar")
    check_arg(sep, "character scalar")
    check_arg(collapse, "NULL character scalar")
    check_arg(frame, "class(environment) l0")
    check_arg(..., "vector len(1)")

    if(...length() == 0) return("")

    sc = sys.call()
    if(identical(sc[[2]], "--help")){
        msg = c(
            "Welcome to dsb help\nUsage: dsb(s) with 's' a character string",
            " ",
            "BASIC usage ------------|",
            "    dsb evaluates anything in '.[]' and inserts it in 's'.",
            '    Ex: if x = "John", then dsb("Hi .[x]!") -> "Hi John!"',
            " ",
            "STRING OPERATIONS ------|",
            "    Each .[] instance supports one or more string operations.",
            "    The syntax is .['arg'op?x] or .['arg'op!x], with:",
            "      - 'arg' a quoted string used as argument,",
            "      - op an operator code,",
            "      - ? or !:",
            "        + ?: evaluates the expression x",
            "        + !: takes x as verbatim",
            "      - x an expression to be evaluated or some verbatim text (no quote needed).",
            '    Ex: dsb(".[\' + \'c?1:3] = 6") -> "1 + 2 + 3 = 6". 1:3 is collapsed (c) with \' + \'.',
            "",
            "    Using ! instead of ? applies the operation to the *verbatim* of the expression.",
            '    Ex: dsb(".[\': => 2\'r!1:3] = 6") -> "123 = 6".',
            "        In the string '1:3', ':' is replaced (r) with '2'.",
            "",
            "    Operations can be chained using comma-separation. The syntax is: .['s1'op1, 's2'op2?x]",
            "    Evaluations are from left to right.",
            '    Ex: dsb(".[\': => 2\'r, \'\'s, \' + \'c!1:3] = 6") -> "1 + 2 + 3 = 6',
            "        1) '1:3'            -> ':' is replaced (r) with '2'  -> '123',",
            "        2) '123'            -> is split (s) with ''          -> c('1', '2', '3')",
            "        3) c('1', '2', '3') -> is collapsed (c) with ' + '   -> '1 + 2 + 3'",
            "",
            "    Nesting works, but only in verbatim components.",
            "    Ex: x = c(\"Doe\", \"Smith\")",
            "        dsb(\"Hi .[' and 'c!John .[x]]\") -> \"Hi John Doe and John Smith\"",
            "",
            "    If nest = TRUE (default), the initial string is implicitly embedded in .[] ",
            "  and considered as verbatim. This means that operations can be applied from the start.",
            "    Ex: dsb(\"', 's!a1, b3, d8\") -> c(\"a1\", \"b3\", \"d8\")",
            "",
            "    Operators have default values, so the quoted argument is optional.",
            '    Ex: dsb("c?1:3") -> "1 2 3". 1:3 is collapsed (c) with \' \', its default.',
            "",
            "OPERATORS --------------|",
            "    Below is the list of operators and their default argument when relevant.",
            "    OPERATOR   DEFAULT                   VALUE ",
            "",
            "    s          ' '                       split, fixed = TRUE",
            "    S          ', *'                     split, perl = TRUE",
            "    x          '[[:alnum:]]+'            extracts the first pattern, perl = TRUE",
            "    X          '[[:alnum:]]+'            extracts all patterns, perl = TRUE",
            "",
            "  Collapse:",
            "    - collapses a vector with paste(x, collapse = 's'),",
            "    - use a double pipe to apply a special collapse to the last values,",
            "    - the syntax is 's1||s2', s2 will be applied to the last 2 values,",
            "    - .[', || and 'c!1:3] -> \"1, 2 and 3\".",
            "    c          ''                        collapse",
            "    C          ', || and '               collapse",
            "",
            "  Replication:",
            "    - use 5* to replicate 5 times,",
            "    - using quotes also works, like in '5'*.",
            "    *          '1'                       replicates n times",
            "    *c         '1'                       replicates n times, then collapses with ''",
            "    **         '1'                       replicates n times, each",
            "    **c        '1'                       replicates n times, each, then collapses with ''",
            "",
            "  Replacement: ",
            "    - the syntax is 'old => new', 'old=>new' or 'old_=>_new',",
            "    - if new is missing, it is considered as the empty string,",
            "    - the default for 'R' is: removing trailing spaces.",
            "    r          '\\n'                      replacement, fixed = TRUE ",
            "    R          '^[ \\t\\r\\n]+|[ \\t\\r\\n]+$' replacement, perl = TRUE ",
            "",
            "  Operators without arguments:",
            "    u                                    puts the first letter of the string to uppercase",
            "    U                                    puts the complete string to uppercase",
            "    l, L                                 puts the complete string to lowercase",
            "    q                                    adds single quotes",
            "    Q                                    adds double quotes",
            "    f                                    applies format(x)",
            "    F                                    applies format(x, justify = 'right')",
            "    w                                    reformat white spaces",
            "    W                                    reformat white spaces and punctuation",
            "    d                                    replaces the content with the empty string",
            "    D                                    removes all elements",
            "    stop                                 removes basic English stop words",
            "",
            "  sprintf formatting:",
            "    - applies a formatting via sprintf,",
            "    - .['.3f'%?x] is equivalent to sprintf('%.3f', x),",
            "    - .['.2f'% ? pi] -> '3.14'.",
            "    %          No default                applies sprintf formatting",
            "",
            "  Keeping a selected number of characters or elements: ",
            "    - the syntax is nK, 'n'K, 'n|s'K, or 'n||s'K",
            "    - with 'n' a number and 's' a string,",
            "    - K keeps the first n elements and drops the rest,",
            "    - k keeps the first n characters of each element,",
            "    - optionaly a string 's' can be appended at the end,",
            "    - but only when the length is greater than n.",
            "    - 'n|s' simply appends 's' while 'n||s' includes the length of 's'",
            "    - in the calculation of the final length.",
            "    - the special markup :rest: and :REST: can be used in 's'.",
            "    - Ex: .['5||..'k!c('hello', 'bonjour')] -> c('hello', 'bon..')",
            "    -     .['5|..'k!c('hello', 'bonjour')]  -> c('hello', 'bonjo..')",
            "    k          NO DEFAULT                keeps the first n characters",
            "    K          NO DEFAULT                keeps the first n elements",
            "    Ko         NO DEFAULT                keeps the first n elements and adds 'others'",
            "                                         ex. usage: 3Ko",
            "",
            "  Insertion and appending operators:",
            "    - the syntax is 's1|s2' to insert, or append, the string s1 first",
            "    - and the string 's2' last.",
            "    - Note that '_|_' can also be used to separate first/last.",
            "    - The uppercase versions insert/append invisibly: ie without changing the vector.",
            "    - for *A, the syntax is 's1|s2 / s3|s4'. If the nbr of elements is:",
            "      * equal to 1: this is equiv. to 's1|s2'A",
            "      * greater then 1: this is equiv. to 's3|s4'A",
            "    a          NO DEFAULT                equiv. to paste0(s1, x, s2)",
            "    A          NO DEFAULT                appends s1/s2 to the beginning/end of the full string",
            "    i          NO DEFAULT                inserts s1/s2 at the beginning/end of the vector",
            "    I          NO DEFAULT                inserts s1/s2 at the beginning/end of the vector",
            "    *A         NO DEFAULT                appends s1/s2 at the beg./end of the full string if len 1",
            "                                         appends s3/s4 if len > 1",
            "",
            "  Adding length dependent verbs/s:",
            "    - verbs can be conjugated WRT the length of the vector,",
            "    - Ex: .[' have'V, C? c('David', 'Steerforth')]  -> 'David and Steerforth have'",
            "          .[' have'V, C? 'Agnes']          -> 'Agnes has'",
            "    - 's' can be added in the same way:",
            "    - Ex: 'The number.[*s_, v, C ? 1:3]' -> 'The numbers are 1, 2 and 3'",
            "          'The number.[*s_, v, C ? 5]'   -> 'The number is 5'",
            "    v          'is '                     adds a conjugated verb at the beginning",
            "    V          ' is'                     adds a conjugated verb at the end",
            "    *s         NO ARGUMENT               adds an 's' at the beginning of the string",
            "    *s_        NO ARGUMENT               adds an 's' followed by a space at the beginning of the string",
            "",
            "CONDITIONS--------------|",
            "    Conditions can be applied with 'if' statements.",
            "    The syntax is 'type comp value'if(true : false), with",
            "    - type: either 'len', 'char', 'fixed' or 'regex'",
            "      + len: number of elements in the vector",
            "      + char: number of characters",
            "      + fixed: fixed pattern",
            "      + regex: regular expression pattern",
            "    - comp: a comparator: ",
            "      + valid for len/char: >, <, >=, <=, !=, ==",
            "      + valid for fixed/regex: !=, ==",
            "    - value: a value for which the comparison is applied.",
            "    - true: operations to be applied if true (can be void)",
            "    - false: operations to be applied if false (can be void)",
            "",
            "    Example: .['char <= 2'if('(|)'a : '[|]'a), ' + 'c ? c(1, 12, 123)]",
            "             ->  '(1) + (12) + [123]'",
            "",
            "  Which operators are valid?",
            "    The 'len' type looks at the length of the vector so that the true/false",
            "    will apply to the full vector. Hence the operations will concern the full vector and",
            "    all operations are valid.",
            "",
            "    The other types apply the true/false operations only to the subset of elements",
            "    which were true or false. Hence to avoid consistency problems, all operators",
            "    that would change the length of the vector/or affect the full vector are disabled.",
            "    Operators not working in that case: A, i, I, K, c, C, *, *c, s, S, v, V.",
            "    Note that the D operator of deletion still works.",
            "",
            "  Special if operators:",
            "    'type comp value'IF(true : false):",
            "      The operator 'IF' works like the regular 'if' except that it also applies any().",
            "      This means that the result of the if statement will be at the level of the vector,",
            "      thus for the types different from 'len', all operators can be applied.",
            "",
            "    *if(true : false):",
            "      The operator *if accepts no argument. The value is true if the length of the vector",
            "      is greater than 1. ",
            "      This is a shortcut to 'len>1'if(true : false).",
            "",
            "SPECIALS ---------------|",
            "    Use '/' first to split the character with commas:",
            '    Ex: dsb("/x1, x2")             -> c("x1", "x2")',
            '        dsb("Hi .[/David, Dora]!") ->  c("Hi David!", "Hi Dora!")',
            "",
            "    In quoted arguments, use backticks to evaluate them from the frame.",
            '    Ex: n = 3 ; dsb("`n`*c!$") -> "$$$". The \'$\' is replicated n times, then collapsed.'
            )

        message(paste(msg, collapse = "\n"))
        return(invisible(NULL))
    }

    res = .dsb(..., frame = frame, sep = sep, vectorize = vectorize, nest = nest, collapse = collapse, check = TRUE)

    class(res) = c("dsb", "character")
    res
}

.dsb0 = function(..., frame = parent.frame(), sep = "", vectorize = FALSE, check = FALSE){
    .dsb(..., frame = frame, nest = FALSE, sep = sep, vectorize = vectorize, check = check)
}

.dsb = function(..., frame = parent.frame(), sep = "", vectorize = FALSE,
                nest = TRUE, collapse = NULL, check = FALSE){

    if(...length() == 0){
        return("")
    } else if(...length() == 1){
        x = as.character(..1)

        if(length(x) > 1){
            stop("dsb can only be applied to character scalars. Problem: the argument is of length ",
                 length(qui), "")
        }

    } else {

        if(check){
            dots = error_sender(list(...), "In dsb, one element of ... could not be evaluated.")
        } else {
            dots = list(...)
        }

        if(any(lengths(dots) > 1)){
            qui = which(lengths(dots) > 1)[1]
            stop("dsb can only be applied to character scalars. Problem: The ", n_th(qui),
                 " elment in ... is of length ", length(dots[[qui]]), ".")
        }

        if(!vectorize){
            # Note: using paste(..1, ..2, sep = sep) explicitly only saves 2us vav do.call
            # not worth it.

            dots$sep = sep
            x = do.call(paste, dots)
        } else {
            # vectorize
            n = length(dots)
            res = vector("list", n)
            for(i in 1:n){
                res[[i]] = .dsb(dots[[i]], nest = nest, frame = frame, check = check)
            }

            return(unlist(res))
        }
    }

    if(is.na(x) || length(x) == 0){
        return(x)
    }

    if(nest){
        x_parsed = list(cpp_dsb_full_string(x))

    } else if(!grepl(".[", x, fixed = TRUE)){
        return(x)

    } else {
        x_parsed = cpp_dsb(x)

    }

    n_x_all = lengths(x_parsed)

    n = length(x_parsed)
    res = character(0)
    i_done = FALSE
    for(i in 1:n){

        if(i_done){
            i_done = FALSE
            next
        }

        xi = x_parsed[[i]]

        if(length(xi) == 1){

            if(i == 1){
                res = xi
            } else {
                res = paste0(res, xi)
            }

        } else {
            operators = xi[[1]]
            xi = xi[[2]]

            if(length(operators) == 0){

                if(nest){
                    # means verbatim => no evaluation
                    if(grepl(".[", xi, fixed = TRUE)){
                        xi = .dsb(xi, frame = frame, nest = FALSE, check = check)
                    }

                } else {
                    # we need to evaluate xi
                    if(check){
                        xi_call = error_sender(str2lang(xi), "The value '", xi,
                                               "' could not be parsed.")
                    } else {
                        xi_call = str2lang(xi)
                    }

                    if(is.character(xi_call)){
                        # if a string literal => it's nesting
                        if(grepl(".[", xi_call, fixed = TRUE)){
                            xi = .dsb(xi_call, frame = frame, nest = FALSE, check = check)
                        }
                    } else {
                        if(check){
                            xi = error_sender(eval(xi_call, frame), "The value '", xi,
                                              "' could not be evaluated.")
                        } else {
                            xi = eval(xi_call, frame)
                        }
                    }
                }

            } else {

                n_op = length(operators)
                verbatim = operators[n_op] %in% c("!", "/")

                if(operators[n_op] == "/"){
                    operators = "', *'S"
                } else {
                    operators = operators[-n_op]
                    # The two separators ? and ! have no default operation
                }

                # If split operator, we concatenate
                concat_nested = length(operators) > 0 && grepl("(s|S)$", operators[[1]])

                if(verbatim && grepl(".[", xi, fixed = TRUE)){
                    xi = .dsb(xi, frame = frame, nest = FALSE, vectorize = concat_nested, check = check)

                } else if(!verbatim){
                    # evaluation
                    if(check){
                        xi_call = error_sender(str2lang(xi), "The value '", xi,
                                               "' could not be parsed.")
                        xi = error_sender(eval(xi_call, frame), "The value '", xi,
                                          "' could not be evaluated.")
                    } else {
                        xi = eval(str2lang(xi), frame)
                    }

                    if(is.function(xi)){
                        stop_up("dsb cannot coerce functions into strings. Problem: '",
                                trimws(x_parsed[[i]][[2]]), "' is a function.")
                    }

                }

                # Now we apply the operators
                for(j in seq_along(operators)){
                    opi = operators[[j]]
                    op_parsed = dsb_char2operator(opi)

                    if(op_parsed$do_eval){
                        if(check){
                            quoted_call = error_sender(str2lang(op_parsed$quoted),
                                                       "In operation '", opi, "', the value '",
                                                       op_parsed$quoted, "' could not be parsed.")

                            quoted = error_sender(eval(quoted_call, frame),
                                                  "In operation '", opi, "', the value '",
                                                  op_parsed$quoted, "' could not be evaluated.")
                        } else {
                            quoted = eval(str2lang(op_parsed$quoted), frame)
                        }

                    } else {
                        quoted = op_parsed$quoted
                    }

                    if(check){
                        xi = error_sender(dsb_operators(xi, quoted, op_parsed$op,
                                                        check = check, frame = frame),
                                          "The operation '", opi, "' failed. Please revise your call.")
                    } else {
                        xi = dsb_operators(xi, quoted, op_parsed$op, check = check, frame = frame)
                    }
                }

                extra = attr(xi, "extra")
                if(length(extra) > 0){

                    if(length(extra$element_add_first) > 0){
                        if(length(extra$element_add_last) > 0){
                            xi = c(extra$element_add_first, xi, extra$element_add_last)
                        } else {
                            xi = c(extra$element_add_first, xi)
                        }
                    } else if(length(extra$element_add_last) > 0){
                        xi = c(xi, extra$element_add_last)
                    }

                    if(length(extra$str_add_first) > 0){
                        if(length(xi) > 0){
                            xi[1] = paste0(extra$str_add_first, xi[1])
                        } else {
                            xi = extra$str_add_first
                        }
                    }

                    if(length(extra$str_add_last) > 0){
                        if(length(xi) > 0){
                            xi[length(xi)] = paste0(xi[length(xi)], extra$str_add_last)
                        } else {
                            xi = extra$str_add_last
                        }
                    }

                    attr(xi, "extra") = NULL
                }
            }

            if(i == 1){
                res = xi
            } else{
                if(i < n && n_x_all[i + 1] == 1){
                    if(vectorize){
                        res = c(res, xi, x_parsed[[i + 1]])
                    } else {
                        res = paste0(res, xi, x_parsed[[i + 1]])
                    }

                    i_done = TRUE
                } else {
                    if(vectorize){
                        res = c(res, xi)
                    } else {
                        res = paste0(res, xi)
                    }
                }
            }
        }
    }

    if(!is.null(collapse) && length(res) > 1){
        res = paste0(res, collapse = collapse)
    }

    return(res)
}


dsb_char2operator = function(x){

    quote = substr(x, 1, 1)

    OPERATORS = c("s", "S", "x", "X", "c", "C", "r", "R",
                  "*", "*c", "**", "**c",
                  "u", "U", "l", "L", "q", "Q", "f", "F", "%",
                  "a", "A", "k", "K", "d", "D", "v", "V",
                  "w", "W", "stop", "i", "I",
                  "*s", "*s_", "*A", "*if", "if", "IF")

    ok = do_eval = FALSE

    if(quote %in% c("'", "\"", "`")){
        pat = paste0("^", quote, "[^", quote, "]*", quote)
        op_abbrev = sub(pat, "", x)
        in_quote = str_trim(x, 1, nchar(op_abbrev) + 1)

        do_eval = quote == "`"

        if(do_eval && substr(in_quote, 1, 1) == "!"){
            # special case: back ticks used as extra quote
            # => enables to use ' and " freely in arguments
            do_eval = FALSE
            in_quote = str_trim(in_quote, 1)
        }

        if(nchar(op_abbrev) == 0){
            stop_up("In dsb, if a quoted value is present, the operators must always be of the form 'value'op, with 'op' an operator. Problem: In '", x, "' the operator is missing.")
        }

    } else if(x %in% OPERATORS){
        # default values
        in_quote = switch(x,
                          s = " ", S = ", *",
                          x = "[[:alnum:]]+", X = "[[:alnum:]]+",
                          c = " ", C = ", || and ",
                          "*" = "1", "*c" = "1",
                          "**" = "1", "**c" = "1",
                          r = "\n", R = "^[ \t\r\n]+|[ \t\r\n]+$",
                          v = "is ", V = " is",
                          "")
        op_abbrev = x

        if(op_abbrev %in% c("%", "k", "K", "a", "A", "if", "IF", "*A")){
            ex = c("%" = ".['.3f'%?pi]",
                   "k" = ".[4k ! longuest word ]", "K" = ".[2K ? 1:5]",
                   "a" = ".['|s'a ! cat]", "A" = ".['|two'A ! one]",
                   "if" = ".['len<3'if(d) ? c('a', 'abc')]",
                   "IF" = ".['char>3'IF(D) ? c('a', 'abc')]",
                   "*A" = ".['The winner/All winners'*A ? c('John', 'Steve')]")
            stop_up("The operator '", op_abbrev,
                    "' has no default value, you must provide values explicitly. Like in ",
                    ex[op_abbrev], " for instance.")
        }

    } else {
        last1 = substr(x, nchar(x), nchar(x))
        if(last1 %in% c("*", "k", "K")){
            in_quote = str_trim(x, -1)

            op_abbrev = last1

            if(op_abbrev == "*" && substr(in_quote, nchar(in_quote), nchar(in_quote)) == "*"){
                op_abbrev = "**"
                in_quote = str_trim(in_quote, -1)
            }

        } else {
            last2 = substr(x, nchar(x) - 1, nchar(x))
            if(last2 %in% c("*c", "Ko", "KO")){
                in_quote = str_trim(x, -2)
                op_abbrev = last2

                if(op_abbrev == "*c" && substr(in_quote, nchar(in_quote), nchar(in_quote)) == "*"){
                    op_abbrev = "**c"
                    in_quote = str_trim(in_quote, -1)
                }

                if(op_abbrev %in% c("Ko", "KO")){
                    # special case
                    text = if(op_abbrev == "Ko") "||:rest: others" else "||:REST: others"
                    in_quote = paste0(in_quote, text)
                    op_abbrev = "K"
                }

            } else {
                first3 = substr(x, 1, 3)

                if(first3 %in% c("*if")){
                    ok = TRUE
                    in_quote = ""
                    op_abbrev = x
                    # OK, dealt with later => special case
                } else {
                    op_abbrev = "problem"
                }
            }
        }
    }

    if(!ok && !op_abbrev %in% OPERATORS){

        first2 = substr(op_abbrev, 1, 2)
        if(!first2 %in% c("if", "IF")){
            msg = c("The operation '", x, "' is not valid. It must be something quoted followed by a valid operator.",
                    "\n  Valid operators: to split: s, S / to replace: r, R  / to collapse: c, C / to extract: x, X",
                    "\n                   to replicate: * / to replicate and collapse with the empty string: *c",
                    "\n                   to upper/lower case: u, U, L / to single/double quote: q, Q",
                    "\n                   to format f, F / to apply sprintf format: %",
                    "\n                   to format whitespaces: w, W / to append: a, A / to insert: i, I",
                    "\n                   to keep: k (#characters), K (#items) / to delete: d, D",
                    "\n                   to remove stopwords: stop / to add a verb: v, V / to add an 's': *s, *s_",
                    "\n                   Conditions: *if(true:false), 'cond'if(true:false), 'cond'IF(true:false)",
                    "\n                   cond: 'regex==expr', 'fixed==expr', 'len==digit', or 'char==digit'.",
                    "\n------------------------------",
                    "\n  type dsb('--help') for more help.",
                    "\n  Example: .[', *'S, 'a => b'r? var] first splits the variable var by commas then replaces every 'a' with a 'b'.")

            message(msg)

            stop_up("In dsb, the operation is not valid, see upper message.")
        }

    }

    res = list(quoted = in_quote, do_eval = do_eval, op = op_abbrev)
    res
}


dsb_operators = function(x, quoted, op, check = FALSE, frame = NULL){

    extra = attr(x, "extra")

    # S, C, R ####

    if(op %in% c("s", "S")){
        # Split is always applied on verbatim stuff => length 1
        if(op == "s"){
            res = unlist(strsplit(x, quoted, fixed = TRUE))
        } else {
            res = unlist(strsplit(x, quoted, perl = TRUE))
        }

        res = res[nchar(res) > 0]

    } else if(op %in% c("c", "C")){
        # collapse

        n_x = length(x)
        if(n_x > 1 && grepl("||", quoted, fixed = TRUE)){
            # This is the "last" operator
            quoted_split = strsplit(quoted, "||", fixed = TRUE)[[1]]
            if(n_x == 2){
                res = paste(x, collapse = quoted_split[[2]])
            } else {
                res = paste(x[-n_x], collapse = quoted_split[[1]])
                res = paste0(res, quoted_split[[2]], x[n_x])
            }
        } else {
            res = paste(x, collapse = quoted)
        }

    } else if(op %in% c("r", "R")){
        new = ""
        if(grepl("=>", quoted, fixed = TRUE)){
            pat = "=>"

            if(grepl(" => ", quoted, fixed = TRUE)){
                pat = " => "
            } else if(grepl("_=>_", quoted, fixed = TRUE)){
                pat = "_=>_"
            }

            quoted_split = strsplit(quoted, pat, fixed = TRUE)[[1]]
            quoted = quoted_split[[1]]
            new = quoted_split[[2]]
        }

        if(op == "r"){
            res = gsub(quoted, new, x, fixed = TRUE)
        } else {
            res = gsub(quoted, new, x, perl = TRUE)
        }

    } else if(op %in% c("*", "**", "*c", "**c")){
        # *, X ####

        if(!is_numeric_in_char(quoted)){
            stop("In dsb: the operator '", op, "' must have numeric arguments, '", quoted, "' is not numeric.")
        }

        if(substr(op, 1, 2) == "**"){
            res = rep(x, each = as.numeric(quoted))
        } else {
            res = rep(x, as.numeric(quoted))
        }

        if(substr(op, nchar(op), nchar(op)) == "c"){
            res = paste(res, collapse = "")
        }

    } else if(op == "x"){
        # extract the first pattern

        x_pat = regexpr(quoted, x, perl = TRUE)

        res = substr(x, x_pat, x_pat - 1 + attr(x_pat, "match.length"))

    } else if(op == "X"){
        # extract all patterns

        x_list = regmatches(x, gregexpr(quoted, x, perl = TRUE))

        res = unlist(x_list)

    } else if(op == "U"){
        # U, L, Q, F, %, W ####

        res = toupper(x)
    } else if(op == "u"){
        # First letter only, if relevant
        res = x
        substr(res, 1, 1) = toupper(substr(x, 1, 1))
    } else if(op %in% c("l", "L")){
        res = tolower(x)
    } else if(op == "q"){
        res = paste0("'", x, "'")
    } else if(op == "Q"){
        res = paste0("\"", x, "\"")
    } else if(op == "f"){
        res = format(x)
        if(is.numeric(x)){
            res = format(trimws(res))
        }
    } else if(op == "F"){
        res = format(x, justify = "right")
    } else if(op == "%"){
        res = sprintf(paste0("%", quoted), x)
    } else if(op %in% c("w", "W")){
        # w: only whitespaces
        # W: whitespaces extended

        res = trimws(x)

        # now the whitespaces
        pat = if(op == "w") "[ \t\r\n]+" else "[[:blank:][:punct:]]+"
        res = gsub(pat, " ", res, perl = TRUE)

    } else if(op %in% c("k", "K")){
        # keep: either the nber of characters (k) or the number of elements (K)

        #
        # Keep ####
        #

        quoted_split = quoted
        pat = c("_||_", "_|_", "||", "|")
        for(p in pat){
            if(grepl(p, quoted, fixed = TRUE)){
                quoted_split = strsplit(quoted, p, fixed = TRUE)[[1]]
                break
            }
        }
        is_included = grepl("||", p, fixed = TRUE)

        if(length(quoted_split) == 1){
            nb = quoted
            add = ""
            is_included = FALSE
        } else {
            nb = quoted_split[[1]]
            add = quoted_split[[2]]
        }

        if(!is_numeric_in_char(nb)){
            stop("In dsb: the operator '", op, "' must first contain a numeric argument, '", quoted, "' does not contain a numeric first.")
        }

        nb = as.numeric(nb)

        if(op == "k"){
            qui = nchar(x) > nb
            res = substr(x, 1, nb)

            if(is_included){
                res[qui] = substr(res[qui], 1, nb - nchar(add))
            }

            if(nchar(add) > 0){
                res[qui] = paste0(res[qui], add)
            }

        } else if(op == "K"){

            if(nb == 0){
                res = character(0)
            } else {
                res = if(nb < length(x)) x[1:nb] else x

                if(any(grepl(":(n|N|rest|REST):", add))){
                    n = length(x)
                    N = ""
                    if(grepl(":N:", add, fixed = TRUE)){
                        N = n_letter(n)
                        add = gsub(":N:", N, add, fixed = TRUE)
                    }

                    add = gsub(":n:", n, add, fixed = TRUE)
                    n_rest = n - nb + is_included
                    if(n_rest > 0){
                        add = gsub(":rest:", n_rest, add, fixed = TRUE)
                        add = gsub(":REST:", n_letter(n_rest), add, fixed = TRUE)
                    }
                }

                if(length(x) > nb){
                    if(is_included){
                        res = res[-nb]
                    }

                    if(nchar(add) > 0){
                        res = c(res, add)
                    }
                }
            }
        }
    } else if(op %in% c("a", "i", "A", "I", "*A")){
        # Appends at the beginning/end of all strings


        # Conditional A
        if(op == "*A"){
            # We transform it into a regular A call
            op = "A"

            pat = c("_/_", " / ", "/")
            quoted_split = quoted
            for(p in pat){
                if(grepl(p, quoted, fixed = TRUE)){
                    quoted_split = strsplit(quoted, p, fixed = TRUE)[[1]]
                    break
                }
            }

            if(length(x) == 1){
                quoted = quoted_split[[1]]
            } else {
                quoted = if(length(quoted_split) > 1) quoted_split[[2]] else ""
            }
        }

        #
        # Insert/Append ####
        #

        pat = c("_|_", "|")
        quoted_split = quoted
        for(p in pat){
            if(grepl(p, quoted, fixed = TRUE)){
                quoted_split = strsplit(quoted, p, fixed = TRUE)[[1]]
                break
            }
        }

        left = quoted_split[[1]]
        right = if(length(quoted_split) > 1) quoted_split[[2]] else ""

        if(nchar(quoted) == 0){
            res = x

        } else if(op == "a"){

            if(length(x) == 0){
                res = x

            } else {

                # We replace the special values
                # :1:, :i:, :a:
                # we allow only one special value
                n_x = length(x)

                pat = c(":1:", ":i:", ":a:", ":I:", ":A:")

                for(i in 1:2){
                    tmp = if(i == 1) left else right
                    any_done = FALSE

                    for(p in pat){
                        if(grepl(p, tmp, fixed = TRUE)){
                            any_done = TRUE
                            tmp_split = strsplit(tmp, p, perl = TRUE)[[1]]
                            if(length(tmp_split) > 1){
                                txt = switch(p,
                                             ":1:" = 1:n_x,
                                             ":i:" = tolower(as.roman(1:n_x)),
                                             ":I:" = as.character(as.roman(1:n_x)),
                                             ":a:" = enum_letter(n_x),
                                             ":A:" = toupper(enum_letter(n_x)))

                                if(length(tmp_split) == 1){
                                    tmp = paste0(tmp_split[[1]], txt)
                                } else  {
                                    tmp_new = paste0(tmp_split[[1]], txt, tmp_split[[2]])
                                    k = 2
                                    while(k + 1 <= length(tmp_split)){
                                        k = k + 1
                                        tmp_new = paste0(tmp_new, txt, tmp_split[[k]])
                                    }
                                    tmp = tmp_new
                                }
                            }

                            break
                        }
                    }

                    if(any_done){
                        if(i == 1){
                            left = tmp
                        } else {
                            right = tmp
                        }
                    }
                }


                if(any(nchar(left) > 0)){
                    if(any(nchar(right) > 0)){
                        res = paste0(left, x, right)
                    } else {
                        res = paste0(left, x)
                    }
                } else {
                    res = paste0(x, right)
                }
            }
        } else if(op == "i"){
            # inserts an ELEMENT at the beginning/end

            if(nchar(left) > 0){
                if(nchar(right) > 0){
                    res = c(left, x, right)
                } else {
                    res = c(left, x)
                }
            } else {
                res = c(x, right)
            }

        } else if(op %in% c("A", "I")){
            # appends/inserts **implicitly** at the beginning of the first/last string

            res = x

            if(any(grepl(":(n|N):", c(left, right)))){
                n = length(x)
                n_letter = ""
                if(any(grepl(":N:", c(left, right), fixed = TRUE))){
                    n_letter = n_letter(n)
                    left = gsub(":N:", n_letter, left)
                    right = gsub(":N:", n_letter, right)
                }

                left = gsub(":n:", n, left)
                right = gsub(":n:", n, right)
            }

            if(is.null(extra)){
                extra = list()
            }

            if(op == "A"){
                if(nchar(left) > 0){
                    extra$str_add_first = paste0(extra$str_add_first, left)
                }

                if(nchar(right) > 0){
                    extra$str_add_last = paste0(extra$str_add_last, right)
                }
            } else if(op == "I"){
                if(nchar(left) > 0){
                    extra$element_add_first = c(extra$element_add_first, left)
                }

                if(nchar(right) > 0){
                    extra$element_add_last = c(extra$element_add_last, right)
                }
            }
        }

        # END: insert/append

    } else if(op == "d"){
        res = rep("", length(x))
    } else if(op == "D"){
        res = character(0)
    } else if(op %in% c("v", "V")){
        # The verb is added implicitly
        # original code from dreamerr

        #
        # *A, D, Verb ####
        #

        res = x

        PLURAL = length(x) > 1

        space_left = space_right = ""
        if(grepl("^ ", quoted)){
            space_left = " "
            quoted = str_trim(quoted, 1)
        }

        if(grepl(" $", quoted)){
            space_right = " "
            quoted = str_trim(quoted, -1)
        }

        past = grepl("\\bpast\\b", quoted)
        if(past){
            quoted = gsub("\\bquoted\\b", "", quoted)
        }

        verb = trimws(quoted)
        if(nchar(verb) == 0){
            verb = "is"
        }

        if(verb %in% c("be", "are")){
            verb = "is"
        } else if(verb == "have"){
            verb = "has"
        } else if(verb == "does"){
            verb = "do"
        } else if(verb %in% c("do not", "does not", "don't", "doesn't")){
            verb = "do not"
        } else if(verb %in% c("is not", "are not", "isn't", "aren't")){
            verb = "is not"
        } else if(verb %in% c("was", "were")){
            verb = "is"
            past = TRUE
        }

        if(past){
            if(verb %in% c("is", "is not", "has", "do", "do not")){
                verb_format = switch(verb, is = ifelse(!PLURAL, "was", "were"), "is not" = ifelse(!PLURAL, "wasn't", "weren't"), has = "had", do = "did", "do not" = "didn't")
            } else {
                verb_format = paste0(verb, "ed")
            }
        } else {
            if(verb %in% c("is", "is not", "has", "do", "do not")){
                verb_format = switch(verb, is = ifelse(!PLURAL, "is", "are"), "is not" = ifelse(!PLURAL, "isn't", "aren't"), has = ifelse(!PLURAL, "has", "have"), do = ifelse(!PLURAL, "does", "do"), "do not" = ifelse(!PLURAL, "doesn't", "don't"))
            } else {
                verb_format = ifelse(PLURAL, verb, paste0(verb, "s"))
            }
        }

        # Adding the verb
        if(is.null(extra)){
            extra = list()
        }

        if(op == "v"){
            extra$str_add_first = paste0(extra$str_add_first, space_left, verb_format, space_right)
        } else if(op == "V"){
            extra$str_add_last = paste0(extra$str_add_last, space_left, verb_format, space_right)
        }

    } else if(op %in% c("*s", "*s_")){

        res = x

        if(is.null(extra)){
            extra = list()
        }

        if(length(x) > 1){
            if(op == "*s_"){
                extra$str_add_first = paste0(extra$str_add_first, "s ")
            } else {
                extra$str_add_first = paste0(extra$str_add_first, "s")
            }

        } else if(op == "*s_"){
            extra$str_add_first = paste0(extra$str_add_first, " ")
        }

    } else if(grepl("^(\\*?if|IF)", op)){
        # The operator is ALWAYS of the form if(yes:no), there must be a parenthesis

        #
        # IF ####
        #

        if(substr(op, 1, 1) == "*"){
            # star if
            quoted = "len > 1"
            if_operators = str_trim(op, 4, 1)
            op = "*if"
        } else {
            if_operators = str_trim(op, 3, 1)
            op = substr(op, 1, 2)
        }

        # We have 4 types:
        # - len
        # - char
        # - fixed
        # - regex

        # left trim
        quoted = sub("^ +", "", quoted)

        type = substr(quoted, 1, 3)
        if(type %in% c("fix", "reg")){
            type = substr(quoted, 1, 5)
        } else if(type == "cha"){
            type = substr(quoted, 1, 4)
        }

        if(!type %in% c("len", "char", "regex", "fixed")){
            stop("In the operator 'if', the argument must be of the form 'type comp value', with 'comp' a comparator (eg ==, >=, etc), and type one of 4 types: 'len', 'char', 'regex' or 'fixed'.",
                 "\nIn the current argument ('", quoted, "') the type is not valid.")
        }

        rest = str_trim(quoted, type)

        is_space = substr(rest, 1, 1) == " "
        if(is_space){
            rest = str_trim(rest, 1)
        }

        comp1 = substr(rest, 1, 1)
        comp2 = substr(rest, 2, 2)

        if(!comp1 %in% c("<", ">", "!", "=") || (comp1 %in% c("!", "=") && comp2 != "=")){
            comp = if(comp1 %in% c("!", "=") && comp2 != "=") paste0(comp1, comp2) else comp1
            stop("In the operator 'if', the argument must be of the form 'type comp value', with 'comp' a comparator (eg ==, >=, etc), and type one of 4 types: 'len', 'char', 'regex' or 'fixed'.",
                 "\nIn the current argument ('", quoted, "') the comparator ('", comp, "') is not valid.")
        }

        if(comp2 == "="){
            comp = paste0(comp1, comp2)
            rest = str_trim(rest, 2)
        } else {
            comp = comp1
            rest = str_trim(rest, 1)
        }

        if(is_space && substr(rest, 1, 1) == " "){
            rest = str_trim(rest, 1)
        }

        value = rest

        # Now we have the 3 elements:
        # - type
        # - comp
        # - value

        if(type %in% c("len", "char")){
            if(!is_numeric_in_char(value)){
                stop("In the operator 'if' (equal to '", quoted, "'), '", type, "' should be compared to a number. Problem: '", value, "' is not a number.")
            }
        }

        is_len = type == "len"

        operators_TF = cpp_dsb_if_extract(if_operators)

        if(isFALSE(operators_TF[[1]])){
            # means error
            stop("Parsing error in the operators of the '", op, "' statement. ",
                 "\nIt must be of the form: ", op, "(yes:no) with yes and no chains of valid operations. Further, 'if' statements cannot be nested.",
                 "\nEx: 'len>3'if(3K, '| and others'A, ', 'c : 'n <= 3'A, D)")
        }

        len_true = FALSE
        if(is_len){
            n_x = length(x)

            len_true = eval(str2lang(paste0(n_x, comp, value)))
            if(len_true){
                x_true = x
                x_false = character(0)
            } else {
                x_true = character(0)
                x_false = x
            }
        } else {
            if(type == "char"){
                comp_call = str2lang(paste0("nchar(x)", comp, value))
            } else {
                comp_txt = paste0("grepl(\"", value, "\", x")
                if(type == "regex"){
                    comp_txt = paste0(comp_txt, ", perl = TRUE)")
                } else {
                    comp_txt = paste0(comp_txt, ", fixed = TRUE)")
                }
                comp_call = str2lang(comp_txt)
            }

            qui = eval(comp_call)

            if(op == "IF"){
                is_len = TRUE
                len_true = any(qui)
                if(len_true){
                    x_true = x
                    x_false = character(0)
                } else {
                    x_true = character(0)
                    x_false = x
                }
            } else {
                x_true = x[qui]
                x_false = x[!qui]
            }
        }

        # All operators that change the length of the element:
        FORBIDDEN = c("c", "C", "*", "*c", "s", "S", "i", "I", "A", "v", "V")

        # if TRUE
        op_true = operators_TF[[1]]

        for(state in c("true", "false")){

            if(state == "true"){
                if(is_len && !len_true) next
                xi = x_true
            } else {
                if(is_len && len_true) next
                xi = x_false
            }

            index = 1 + (state == "false")
            op_all = operators_TF[[index]]

            for(i in seq_along(op_all)){
                opi = op_all[i]

                op_parsed = dsb_char2operator(opi)

                if(!is_len && op_parsed$op %in% FORBIDDEN){
                    stop_up("In the 'if' statement, you cannot use operators that would change the length of the vector, hence '", op_parsed$op, "' is forbidden. The IF statement (which does like 'if(any(cond))'), or if applied to the length, would be OK since they apply changes to the full vector.")
                }

                if(op_parsed$do_eval){
                    if(check){
                        quoted_call = error_sender(up = 1, str2lang(op_parsed$quoted),
                                                   "In operation '", opi, "', the value '",
                                                   op_parsed$quoted, "' could not be parsed.")

                        quoted = error_sender(up = 1, eval(quoted_call, frame),
                                              "In operation '", opi, "', the value '",
                                              op_parsed$quoted, "' could not be evaluated.")
                    } else {
                        quoted = eval(str2lang(op_parsed$quoted), frame)
                    }

                } else {
                    quoted = op_parsed$quoted
                }

                if(check){
                    xi = error_sender(up = 1, dsb_operators(xi, quoted, op_parsed$op),
                                      "The operation '", opi, "' failed. Please revise your call.")
                } else {
                    xi = dsb_operators(xi, quoted, op_parsed$op)
                }
            }

            if(state == "true"){
                x_true = xi
            } else {
                x_false = xi
            }
        }

        if(is_len){
            res = if(len_true) x_true else x_false

            # we can add extra stuff, if needed

            extra_if = attr(res, "extra")
            if(length(extra_if) > 0){
                if(is.null(extra)){
                    extra = list()
                }

                if(length(extra_if$element_add_first) > 0){
                    extra$element_add_first = c(extra$element_add_first, extra_if$element_add_first)
                }

                if(length(extra_if$element_add_last) > 0){
                    extra$element_add_last = c(extra$element_add_last, extra_if$element_add_last)
                }

                if(length(extra_if$str_add_first) > 0){
                    extra$str_add_first = paste0(extra$str_add_first, extra_if$str_add_first)
                }

                if(length(extra_if$str_add_last) > 0){
                    extra$str_add_last = paste0(extra$str_add_last, extra_if$str_add_last)
                }
            }

        } else {
            # extra is not touched

            res = x

            if(length(x_true) == 0 && length(x_false) == 0){
                res = character(0)
            } else if(length(x_true) == 0){
                res[!qui] = x_false
                res = res[!qui]
            } else if(length(x_false) == 0){
                res[qui] = x_true
                res = res[qui]
            } else {
                res[qui] = x_true
                res[!qui] = x_false
            }
        }

    } else if(op == "stop"){
        # stop ####

        # current limitation: does not work for quoted words
        #                     but quoted words are not stopwords usually

        # Snowball stopwords
        # These come from http://snowballstem.org/algorithms/english/stop.txt
        stopwords = c("i", "me", "my", "myself", "we", "our", "ours", "ourselves",
                      "you", "your", "yours", "yourself", "yourselves", "he", "him",
                      "his", "himself", "she", "her", "hers", "herself", "it", "its",
                      "itself", "they", "them", "their", "theirs", "themselves", "what",
                      "which", "who", "whom", "this", "that", "these", "those", "am", "is",
                      "are", "was", "were", "be", "been", "being", "have", "has", "had",
                      "having", "do", "does", "did", "doing", "would", "should",
                      "could", "ought", "i'm", "you're", "he's", "she's", "it's",
                      "we're", "they're", "i've", "you've", "we've", "they've",
                      "i'd", "you'd", "he'd", "she'd", "we'd", "they'd", "i'll",
                      "you'll", "he'll", "she'll", "we'll", "they'll", "isn't",
                      "aren't", "wasn't", "weren't", "hasn't", "haven't", "hadn't",
                      "doesn't", "don't", "didn't", "won't", "wouldn't", "shan't",
                      "shouldn't", "can't", "cannot", "couldn't", "mustn't", "let's",
                      "that's", "who's", "what's", "here's", "there's", "when's", "where's",
                      "why's", "how's", "a", "an", "the", "and", "but", "if", "or",
                      "because", "as", "until", "while", "of", "at", "by", "for",
                      "with", "about", "against", "between", "into", "through",
                      "during", "before", "after", "above", "below", "to", "from",
                      "up", "down", "in", "out", "on", "off", "over", "under",
                      "again", "further", "then", "once", "here", "there", "when",
                      "where", "why", "how", "all", "any", "both", "each", "few",
                      "more", "most", "other", "some", "such", "no", "nor", "not",
                      "only", "own", "same", "so", "than", "too", "very")

        n = length(x)
        x_split = strsplit(x, "(?<=[[:alnum:]])(?=[^[:alnum:]'])|(?<=[^[:alnum:]'])(?=[[:alnum:]])", perl = TRUE)
        x_len = lengths(x_split)
        x_vec = unlist(x_split)

        id = rep(1:n, x_len)

        # Lowering is costly, checking costs a bit less
        if(any(grepl("[[:upper:]]", x))){
            qui_drop = which(tolower(x_vec) %in% stopwords)
        } else {
            qui_drop = which(x_vec %in% stopwords)
        }

        x_vec = x_vec[-qui_drop]
        id = id[-qui_drop]

        res = cpp_paste_conditional(x_vec, id, n)


    } else {
        stop("In dsb: the operator '", op, "' is not recognized. Internal error: this problem should have been spotted beforehand.")
    }

    if(length(extra) > 0){
        attr(res, "extra") = extra
    }

    return(res)
}




n_letter = function(n){
    num2char = c("zero", "one", "two", "three", "four", "five", "six", "seven",
                 "eight", "nine", "ten", "eleven", "twelve", "thirteen",
                 "fourteen", "fifteen", "sixteen", "seventeen", "eighteen",
                 "nineteen")
    if(n < 20){
        n_letter = num2char[n + 1]
    } else if(n < 100){
        tens = n %/% 10
        digit = n %% 10
        tens_letter = c("twenty", "thirty", "forty", "fifty",
                        "sixty", "seventy", "eighty", "ninety")

        num2char = paste0("-", num2char)
        num2char[1] = ""

        n_letter = paste0(tens_letter[tens - 1], num2char[digit + 1])
    } else {
        n_letter = n
    }

    n_letter
}

enum_letter = function(n){
    # returns only lowercase
    # oddity: there is no powers of z
    # => because otherwise the algorithm would have been too complex
    # and we don't care tbh


    if(n < 27){
        return(letters[1:n])
    }

    n_26 = log(n, 26)
    n_26 = floor(n_26) + (n_26 %% 1 == 0)

    rest = (1:n) - 1
    res = vector("list", n_26 + 1)
    i = 1
    for(p in n_26:0){
        num = 26**p

        if(p == 0){
            res[[i]] = letters[(rest %/% num) + 1]
        } else {
            res[[i]] = c("", letters)[(rest %/% num) + 1]
        }

        rest = rest %% num
        i = i + 1
    }

    res = do.call(base::paste0, res)

    res
}



