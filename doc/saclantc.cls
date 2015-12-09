%% 
%% This is file `saclantc.cls', generated on <1995/11/13> 
%% with the docstrip utility (2.2i).
%% 
%% The original source file was `saclantc.dtx'.
%% 
%% 
%% IMPORTANT NOTICE:
%% You are not allowed to distribute this file.
%% For distribution of the original source see
%% the copyright notice in the file saclantc.dtx.
%% 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{saclantc}
             [1995/10/23 v0.2 updatad 22 May 1997
             `saclantc' LaTeX document class]
\newcommand{\SC@class}{}
\newif\if@SC@camera
\@SC@camerafalse
\newif\if@restonecol
\newif\if@openright
\newif\if@openbib       \@openbibfalse
\newif\if@outsideauthor \@outsideauthorfalse
\newif\if@mainmatter    \@mainmattertrue
%\DeclareOption{nu}{\renewcommand{\SC@class}{NATO UNCLASSIFIED}}
\DeclareOption{nu}{\renewcommand{\SC@class}{}}
\DeclareOption{nr}{\renewcommand{\SC@class}{NATO RESTRICTED}}
\DeclareOption{nc}{\renewcommand{\SC@class}{NATO CONFIDENTIAL}}
\DeclareOption{ns}{\renewcommand{\SC@class}{NATO SECRET}}
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\@SC@camerafalse
  \setlength\overfullrule{5\p@}\@twosidefalse \@mparswitchfalse}
\DeclareOption{final}{\@SC@cameratrue
  \setlength\overfullrule{\z@}\@twosidetrue  \@mparswitchtrue}
\DeclareOption{openright}{\@openrighttrue}
\DeclareOption{openany}{\@openrightfalse}
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{\@openbibtrue}
\ExecuteOptions{nu,a4paper,11pt,oneside,onecolumn,draft}
\ProcessOptions
\input{size11.clo}
\def\magstep
        #1{\ifcase#1 1000\or1200\or1440\or1728
        \or2074\or2488\or2986\fi\relax}

\def\@magscale#1{ scaled \magstep #1}
\def\@halfmag{ scaled \magstephalf}

\font\ten@ss@I=cmss10 \@magscale1        % for titlepage
\font\ten@ss@II=cmss10 \@magscale2       % for cover and title pages
\font\ten@ss@III=cmss10 \@magscale3      % for coverpage

\font\ten@ssi@II=cmssi10 \@magscale2     % for cover and title pages
\font\ten@ssi@III=cmssi10 \@magscale3    % for coverpage
\font\ten@ssi@IV=cmssi12 \@magscale4     % for headings
\font\ten@ssi@VI=cmssi12 \@magscale6     % for headings

\font\ten@ssbf@=cmssbx10 \@halfmag       % for Distribution string
\font\ten@ssbf@II=cmssbx10 \@magscale2   % for coverpage
\font\ten@ssbf@IV=cmssbx10 \@magscale4   % for cover classification
\font\ten@ssbf@V=cmbx10 \@magscale5      % for vugraph titles

\font\eight@ssq@=cmssq8 \@halfmag        % for tagline classification
\font\eight@ssqi@=cmssqi8 \@halfmag      % for tagline tags
\font\eight@ssq@I=cmssq8 \@magscale1     %
\font\eight@ssqi@I=cmssqi8 \@magscale1   % for subheadings

\font\twelve@rm@=cmr12 \@halfmag         % for \SC@stamp cap ess
\font\eight@rm@=cmr8 \@halfmag           % for \SC@stamp small caps
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand{\baselinestretch}{}
\setlength\parskip{12\p@ \@plus \p@}
\setlength\parindent{\z@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setlength\headsep   {35\p@}
\setlength\topskip   {11\p@}
\setlength\footskip{35\p@}
\setlength\maxdepth{.5\topskip}
\setlength\@maxdepth\maxdepth
\setlength\@tempdima{\paperwidth}
\addtolength\@tempdima{-2in}
\setlength\@tempdimb{360\p@}
\if@twocolumn
  \ifdim\@tempdima>2\@tempdimb\relax
    \setlength\textwidth{2\@tempdimb}
  \else
    \setlength\textwidth{\@tempdima}
  \fi
\else
  \ifdim\@tempdima>\@tempdimb\relax
    \setlength\textwidth{\@tempdimb}
  \else
    \setlength\textwidth{\@tempdima}
  \fi
\fi
\@settopoint\textwidth
\setlength\@tempdima{\paperheight}
\addtolength\@tempdima{-2in}
\addtolength\@tempdima{-1.5in}
\divide\@tempdima\baselineskip
\@tempcnta=\@tempdima
\setlength\textheight{\@tempcnta\baselineskip}
\addtolength\textheight{\topskip}
\setlength\textwidth{14cm}
\setlength\textheight{21cm}
\if@twoside
  \setlength\@tempdima        {\paperwidth}
  \addtolength\@tempdima      {-\textwidth}
  \setlength\oddsidemargin    {.4\@tempdima}
  \addtolength\oddsidemargin  {-1in}
  \setlength\marginparwidth   {.6\@tempdima}
  \addtolength\marginparwidth {-0.4in}
\else
  \setlength\@tempdima        {\paperwidth}
  \addtolength\@tempdima      {-\textwidth}
  \setlength\oddsidemargin    {.5\@tempdima}
  \addtolength\oddsidemargin  {-1in}
  \setlength\marginparwidth   {.5\@tempdima}
  \addtolength\marginparwidth {-.4in}
\fi
\ifdim \marginparwidth >2in
   \setlength\marginparwidth{2in}
\fi
\@settopoint\oddsidemargin
\@settopoint\marginparwidth
\setlength\evensidemargin  {\paperwidth}
\addtolength\evensidemargin{-2in}
\addtolength\evensidemargin{-\textwidth}
\addtolength\evensidemargin{-\oddsidemargin}
\@settopoint\evensidemargin
\if@twocolumn
 \setlength\marginparsep {10\p@}
\else
  \setlength\marginparsep{10\p@}
\fi
\setlength\marginparpush{5\p@}
\setlength\topmargin{\paperheight}
\addtolength\topmargin{-2in}
\addtolength\topmargin{-\headheight}
\addtolength\topmargin{-\headsep}
\addtolength\topmargin{-\textheight}
\addtolength\topmargin{-\footskip}     % this might be wrong!
\addtolength\topmargin{-.5\topmargin}
\@settopoint\topmargin
\voffset=-1in \hoffset=-1in        % let's zero the origin
\setlength{\topmargin}{3cm}
\setlength{\oddsidemargin}{3.5cm}
\setlength{\evensidemargin}{3.5cm}
\setlength\footnotesep{7.7\p@}
\setlength{\skip\footins}{10\p@ \@plus 4\p@ \@minus 2\p@}
\setcounter{topnumber}{2}
\renewcommand{\topfraction}{.7}
\setcounter{bottomnumber}{1}
\renewcommand{\bottomfraction}{.3}
\setcounter{totalnumber}{3}
\renewcommand{\textfraction}{.2}
\renewcommand{\floatpagefraction}{.5}
\setcounter{dbltopnumber}{2}
\renewcommand{\dbltopfraction}{.7}
\renewcommand{\dblfloatpagefraction}{.5}
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\@fptop{0\p@ \@plus 1fil}
\setlength\@fpsep{8\p@ \@plus 2fil}
\setlength\@fpbot{0\p@ \@plus 1fil}
\setlength\@dblfptop{0\p@ \@plus 1fil}
\setlength\@dblfpsep{8\p@ \@plus 2fil}
\setlength\@dblfpbot{0\p@ \@plus 1fil}
\if@twoside
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@evenhead{\thepage\hfil\slshape\leftmark}%
      \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {\uppercase{%
        \ifnum \c@secnumdepth >\z@
          \thesection.\hskip 1em\relax
        \fi
        ##1}}{}}%
    \def\subsectionmark##1{%
      \markright {%
        \ifnum \c@secnumdepth >\@ne
          \thesubsection\hskip 1em\relax
          \fi
          ##1}}}
\else
  \def\ps@headings{%
    \let\@oddfoot\@empty
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markright {\uppercase{%
        \ifnum \c@secnumdepth >\m@ne
          \thesection\hskip 1em\relax
        \fi
        ##1}}}}
\fi
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
    }
\newbox\SC@stamp
\def\set@SC@stamp
        {%
        \setbox\SC@stamp\vbox{\baselineskip12\p@
        \setbox\z@\hbox{\twelve@rm@ S} \dimen@\wd\z@
        \setbox\tw@\hbox{\eight@rm@ ACLANTCEN
        \SC@serialno\unskip}\advance\dimen@ by\wd\tw@
        \hbox to\dimen@{\kern\wd\z@\hrulefill}\kern-4\p@
        \hbox to\dimen@{\box\z@\box\tw@}%
        \kern-8.6\p@\hbox to\dimen@{\hrulefill}\kern-3\p@}%
        }
\AtBeginDocument{\set@SC@stamp}
\def\SC@class@set{{\sffamily\bfseries\SC@class}}

\newcommand{\ps@class}{%
  \renewcommand{\@oddhead}{\copy\SC@stamp\hfil\SC@class@set}%
  \renewcommand{\@evenhead}{\SC@class@set\hfil\copy\SC@stamp}%
  \renewcommand{\@oddfoot}%
      {\ifnum\c@page=\m@ne\hfil\SC@class@set\else
      \hss\textrm{-- \thepage\ --}\hss\llap\SC@class@set\fi}%
  \renewcommand{\@evenfoot}%
      {\rlap\SC@class@set\hss\textrm{-- \thepage\ --}\hss}%
}
\newcommand*{\division}[1]{\gdef\SC@division{#1}}
\newcommand*{\project}[1]{\gdef\SC@project{#1}}
\newcommand*{\serialno}[1]{\gdef\SC@serialno{#1}}
\newcommand*{\director}[1]{\gdef\@director{#1}}
\newcommand*{\authoroutside}[1]{\@outsideauthortrue
                               \gdef\@authoroutside{#1}}
\newcommand*{\outsideinfo}[1]{\gdef\@outsideinfo{#1}}
\def\NATO{North Atlantic Treaty Organization}
\def\SACLANTCEN{SACLANT Undersea Research Centre}
\hyphenation{SAC-LANT-CEN}

\def\div@chief
        {\ifx\doc@division\@SRD
        \else\ifx\doc@division\@ERD
        \fi\fi}
\def\@division
        {\ifx\doc@division\@SRD
        Systems Research Division\else\ifx\doc@division\@ERD
        Environmental Research Division\fi\fi}

\def\released@by
        {\ifx\doc@type\CEN@REPT
                \let\@release\empty
                \def\@signed{\@director\break Director}
        \else\ifx\doc@type\CEN@MEMO
                \def\@release{Issued by:\break\@division}
                \def\@signed{\div@chief\break Division Chief}
        \else\ifx\doc@type\CEN@SPEC
                \def\@release{Issued by:\break\@division}
                \def\@signed{\div@chief\break Division Chief}
        \else\ifx\doc@type\PRO@WORK
                \let\@release\empty
                \def\@signed{\@director\break Director}\fi\fi\fi\fi}

                \let\@release\empty
                \def\@signed{Director}
\def\SC@releasedby{Director}
\def\SC@mumble
 {}
  % peter       {The content of this document pertains to work performed under
 % peter        Project \SC@project\ of the SACLANTCEN Programme of Work.
  % peter       The document has been approved for release by The Director,
   % peter      SACLANTCEN.}
\def\SC@nohyphens
        {\pretolerance=\@M \tolerance=\@M
        \hyphenpenalty=\@M \exhyphenpenalty=\@M}

\def\SC@center
        {\leftskip=\z@ plus4em \rightskip=\leftskip
        \parfillskip\z@ \spaceskip=0.3333em \xspaceskip=0.5em}

\def\dlap#1{\vbox to\z@{#1\vss}}
\def\ulap#1{\vbox to\z@{\vss#1}}
\def\udlap#1{\vbox to\z@{\vss#1\vss}}
\def\lrlap#1{\hbox to\z@{\hss#1\hss}}

\newcommand{\maketitle}{\begin{titlepage}%
        \vbox to\textheight
        {\SC@nohyphens
%        \centerline{\dlap{\hsize7.5cm
        \centerline{\dlap{\hsize10cm
                \vskip7.5cm
                {\center{\ulap{\SC@center\ten@ss@II\baselineskip16\p@\@title}
                 }}
                \vskip1cm
                \dlap{{\SC@center\ten@ss@I\baselineskip14\p@\@author
                         \vskip1cm
                         \nointerlineskip
                        \centerline{\vrule width2.5cm height.4\p@}
                        \nointerlineskip}
                        \vskip.5cm\centerline
                         {\dlap{\hsize6.35cm\parfillskip\z@\SC@mumble}}
                        }
                }}
 % peter       \vfill
 % peter       \released@by
 % peter       \centerline
  % peter              {\ulap{\SC@center\eight@ssq@\@release\vskip3cm\@signed}}
        }%
% peter  \clearpage
%  peter \vspace*{\fill}
% peter  \begin{center}
%  peter   intentionally blank page
%  peter \end{center}
%  peter \vfill
  \end{titlepage}%
  \let\maketitle\relax
}
\setcounter{secnumdepth}{3}
\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand{\thepart}         {\Roman{part}}
\renewcommand{\thesection}      {\arabic{section}}
\renewcommand{\thesubsection}   {\thesection.\arabic{subsection}}
\renewcommand{\thesubsubsection}{\thesubsection.\arabic{subsubsection}}
\renewcommand{\theparagraph}    {\thesubsubsection.\arabic{paragraph}}
\renewcommand{\thesubparagraph} {\theparagraph.\arabic{subparagraph}}
\newcommand{\@sectannx}{\sectionname}
\newcommand{\frontmatter}{%
            \@mainmatterfalse\pagenumbering{roman}}
\newcommand{\mainmatter}{\cleardoublepage
       \@mainmattertrue\pagenumbering{arabic}}
\newcommand{\backmatter}{\if@openright\cleardoublepage\else\clearpage\fi
      \@mainmatterfalse}
\newenvironment{prelims}{%
    \frontmatter\small\maketitle\titlepage
 % peter   \cleardoublepage
}
    {\par\vfil\null\endtitlepage\mainmatter}
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{class}%
      }%
    {\if@restonecol\twocolumn \else \cleardoublepage \fi
    }
\def\SC@title@and@author
        {\vglue1cm
        \vtop{\hsize6cm\rightskip\z@ plus2cm minus2cm
        \parskip\z@ \let\break\newline
        \SC@nohyphens \bfseries\@title \medskip
        \hrule width1.5cm \medskip
        \mdseries\@author \bigskip}}
\def\SC@title@and@authoroutside
        {\vglue1cm
        \vtop{\hsize6cm\rightskip\z@ plus2cm minus2cm
        \parskip\z@ \let\break\newline
        \SC@nohyphens \bfseries\@title \medskip
        \hrule width1.5cm \medskip
        \mdseries\@authoroutside \bigskip}}
\newenvironment{executive}{%
    \small 
    \if@outsideauthor
       \SC@title@and@authoroutside
    \else
       \SC@title@and@author
    \fi
    \leftskip=1cm \rightskip=\leftskip
      {\bfseries \executivename:\quad}}
   {
			\vfill 
  \if@outsideauthor
    {\small \@outsideinfo}
  \fi
% peter  \clearpage
% peter   \vspace*{\fill}
% peter   \begin{center}
% peter     intentionally blank page
% peter   \end{center}
% peter   \vfill
% peter   \cleardoublepage
}
\newenvironment{abstract}{%
    \small \SC@title@and@author
    \leftskip=1cm \rightskip=\leftskip
      {\bfseries \abstractname:\quad}}
   {\clearpage}
\newenvironment{keywords}{%
    \small
    \leftskip=1cm \rightskip=\leftskip
      {\bfseries \keywordsname:\quad}}
   {\clearpage}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART HEADING
\newcommand{\part}{\cleardoublepage
                 \thispagestyle{class}%
                 \if@twocolumn
                     \onecolumn
                     \@tempswatrue
                   \else
                     \@tempswafalse
                 \fi
                 \hbox{}\vfil
                 \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    \markboth{}{}
    {\centering
     \interlinepenalty \@M
     \reset@font
     \ifnum \c@secnumdepth >-2\relax
       \huge\bfseries \partname~\thepart
       \par
       \vskip 20\p@
     \fi
     \Huge \bfseries #2\par}%
    \@endpart}
\def\@spart#1{%
    {\centering
     \interlinepenalty \@M
     \reset@font
     \Huge \bfseries #1\par}%
    \@endpart}
\def\@endpart{\vfil\newpage
              \if@twoside
                \hbox{}%
                \thispagestyle{class}%
                \newpage
              \fi
              \if@tempswa
                \twocolumn
              \fi}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION HEADING
\newcommand{\section}{\if@openright\cleardoublepage\else\clearpage\fi
                    \thispagestyle{class}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@section\@ssection}
\def\@section[#1]#2{\ifnum \c@secnumdepth >\m@ne
                         \refstepcounter{section}%
                         \typeout{\@sectannx\space\thesection.}%
                         \addcontentsline{toc}{section}%
                                   {\protect\numberline{\thesection}#1}%
                    \else
                      \addcontentsline{toc}{section}{#1}
                    \fi
                    \sectionmark{#1}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makesectionhead{#2}]%
                    \else
                      \@makesectionhead{#2}%
                      \@afterheading
                    \fi}
\def\SC@raggedleft
       {\ifodd\c@page\leftskip.65\textwidth
               plus.35\textwidth minus1.0\textwidth
       \else\rightskip.65\textwidth
               plus.35\textwidth minus1.0\textwidth
       \fi
       \parskip\z@\parfillskip\z@\hyphenpenalty\@M}

\global\let\@parskip\parskip
\def\SC@nobreak
       #1{\parskip\z@\nobreak#1\gdef\par{\endgraf\parskip\@parskip
       \global\let\par\endgraf}}

\def\@makesectionhead#1{%
  \vspace*{.05\textheight}%
  \reset@font
  {\parindent \z@
        \vtop{\baselineskip16\p@\SC@raggedleft
        \hbadness\@M
        \ifnum \c@secnumdepth >\m@ne
          {\ten@ssi@VI \thesection}%
        \fi
        \vadjust{\vskip6\p@}\linebreak
        \Large\sffamily #1\linebreak
        \raise6\p@\hbox to\hsize{\hrulefill}}\par
        \SC@nobreak\smallskip
  }}
\def\@ssection#1{\if@twocolumn
                   \@topnewpage[\@makessectionhead{#1}]%
                 \else
                   \@makessectionhead{#1}%
                   \@afterheading
                 \fi}
\def\@makessectionhead#1{%
  \vspace*{.05\textheight}%
  \reset@font
  {\parindent \z@
        \vtop{\baselineskip16\p@\SC@raggedleft
        \hbadness\@M
        \Large\sffamily #1\linebreak
        \raise6\p@\hbox to\hsize{\hrulefill}}\par
        \SC@nobreak\smallskip
  }}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANNEX HEADING
\newcommand{\annex}{\par
  \if@mainmatter\backmatter
    \setcounter{section}{0}%
    \renewcommand{\@sectannx}{\annexname}%
    \renewcommand{\thesection}{\Alph{section}}%
  \fi
  \annexsection}
\newcommand{\annexsection}{\if@openright\cleardoublepage\else\clearpage\fi
                    \thispagestyle{class}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \@dblarg{\@annexsection}}
%                    \secdef\@section\@ssection}
\def\@annexsection[#1]#2{\ifnum \c@secnumdepth >\m@ne
                         \refstepcounter{section}%
                         \typeout{\@sectannx\space\thesection.}%
                         \addcontentsline{toc}{annex}%
                                   {\protect\numberline{Annex~\thesection~-~}#1}%
                    \else
                      \addcontentsline{toc}{annex}{#1}
                    \fi
                    \sectionmark{#1}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makeannexhead{#2}]%
                    \else
                      \@makeannexhead{#2}%
                      \@afterheading
                    \fi}
\def\@makeannexhead#1{%
  \vspace*{.05\textheight}%
  \reset@font
  {\parindent \z@
        \vtop{\baselineskip16\p@\SC@raggedleft
        \hbadness\@M
        \ifnum \c@secnumdepth >\m@ne
          {\Large\sffamily Annex \thesection}%
        \fi
        \vadjust{\vskip6\p@}\linebreak
        \Large\sffamily #1\linebreak
        \raise6\p@\hbox to\hsize{\hrulefill}}\par
        \SC@nobreak\smallskip
  }}

%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\subsection}{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\reset@font\large\eight@ssqi@I}}
\newcommand{\subsubsection}{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\reset@font\normalsize\eight@ssqi@I}}
\newcommand{\paragraph}{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\reset@font\normalsize}}
\newcommand{\subparagraph}{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\reset@font\normalsize}}
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength\leftmargin    {\leftmargini}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\setlength\partopsep{3\p@ \@plus 1\p@ \@minus 1\p@}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand{\theenumi}{\arabic{enumi}}
\renewcommand{\theenumii}{\alph{enumii}}
\renewcommand{\theenumiii}{\roman{enumiii}}
\renewcommand{\theenumiv}{\Alph{enumiv}}
\newcommand{\labelenumi}{\theenumi.}
\newcommand{\labelenumii}{(\theenumii)}
\newcommand{\labelenumiii}{\theenumiii.}
\newcommand{\labelenumiv}{\theenumiv.}
\renewcommand{\p@enumii}{\theenumi}
\renewcommand{\p@enumiii}{\theenumi(\theenumii)}
\renewcommand{\p@enumiv}{\p@enumiii\theenumiii}
\newcommand{\labelitemi}{$\m@th\bullet$}
\newcommand{\labelitemii}{\normalfont\bfseries --}
\newcommand{\labelitemiii}{$\m@th\ast$}
\newcommand{\labelitemiv}{$\m@th\cdot$}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand{\descriptionlabel}[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
\newenvironment{verse}
               {\let\\=\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item[]}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item[]}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item[]}
               {\endlist}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand{\theequation}{\arabic{equation}}
\newcounter{figure}
\renewcommand{\thefigure}{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\bf #1}\hspace{0.6em} {\it #2}}%
  \ifdim \wd\@tempboxa >\hsize
    {\bf #1}\hspace{0.6em} {\it #2}\par
  \else
    \global \@minipagefalse
    \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*{\cal}{\@fontswitch{\relax}{\mathcal}}
\DeclareRobustCommand*{\mit}{\@fontswitch{\relax}{\mathnormal}}
\newcommand{\@pnumwidth}{1.55em}
\newcommand{\@tocrmarg} {2.55em}
\newcommand{\@dotsep}{4.5}
\setcounter{tocdepth}{2}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE OF CONTENTS
\newcommand{\tableofcontents}{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \tocsection{\contentsname
        \@mkboth{\uppercase{\contentsname}}{\uppercase{\contentsname}}}%
    \small
    \leftskip=1cm \rightskip=\leftskip
    \parskip\z@
    \@starttoc{toc}%
    \if@restonecol\twocolumn\fi
    }
\newcommand{\tocsection}{\if@openright\cleardoublepage\else\clearpage\fi
                    \thispagestyle{class}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \@tocssection}
\def\@tocssection#1{\if@twocolumn
                   \@topnewpage[\@maketochead{#1}]%
                 \else
                   \@maketochead{#1}%
                   \@afterheading
                 \fi}
\def\@maketochead#1{%
  \vspace*{.05\textheight}%
  \reset@font
  {\parindent \z@
        \vtop{\baselineskip16\p@\SC@raggedleft
        \hbadness\@M
        \Large\sffamily #1\linebreak}
        \SC@nobreak\smallskip
  }}
\newcommand*{\l@part}[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{\@secpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength\@tempdima{3em}%
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
       \nobreak
    \endgroup
  \fi}
\newcommand*{\l@section}[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty{\@secpenalty}%
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parskip \z@
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode
      \advance\leftskip\@tempdima
      \hskip -\leftskip
%      #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
% Want dotted TOC line for sections (DAA, 2/97)
    {#1}\nobreak\leaders\hbox{$\m@th \mkern \@dotsep mu.\mkern \@dotsep
       mu$}\hfill \nobreak
           \hbox to\@pnumwidth{%
             \hfil\reset@font\rmfamily \normalcolor #2}\par
    \endgroup
  \fi}

\newcommand*{\l@references}[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty{\@secpenalty}%
    \addvspace{1.5em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parskip \z@
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode
      \advance\leftskip\@tempdima
      \hskip -\leftskip
%      #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
% Want dotted TOC line for sections (DAA, 2/97)
    {#1}\nobreak\leaders\hbox{$\m@th \mkern \@dotsep mu.\mkern \@dotsep
       mu$}\hfill \nobreak
           \hbox to\@pnumwidth{%
             \hfil\reset@font\rmfamily \normalcolor #2}\par\vspace{1.5em}
    \endgroup
  \fi}

\newcommand*{\l@annex}        {\@dottedtocline{1}{0em}{5em}}
\newcommand*{\l@subsection}   {\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*{\l@subsubsection}{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*{\l@paragraph}    {\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*{\l@subparagraph} {\@dottedtocline{5}{10em}{5em}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST OF FIGURES
\newcommand{\listoffigures}{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \section*{\listfigurename
      \@mkboth{\uppercase{\listfigurename}}%
              {\uppercase{\listfigurename}}}%
    \@starttoc{lof}%
    \if@restonecol\twocolumn\fi
    }
\newcommand*{\l@figure}{\@dottedtocline{1}{1.5em}{2.3em}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST OF TABLES
\newcommand{\listoftables}{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \section*{\listtablename
      \@mkboth{\uppercase{\listtablename}}{\uppercase{\listtablename}}}%
    \@starttoc{lot}%
    \if@restonecol\twocolumn\fi
    }
\let\l@table\l@figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BIBLIOGRAPHY
\newdimen\bibindent
\bibindent=1.5em
\newenvironment{thebibliography}[1]
     {\section*{\refname
        \@mkboth{\uppercase{\refname}}{\uppercase{\refname}}}%
      \addcontentsline{toc}{references}{\refname}  		    % Want references in TOC (DAA, 2/97)
      \list{\@biblabel{\arabic{enumiv}}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \if@openbib
              \advance\leftmargin\bibindent
              \itemindent -\bibindent
              \listparindent \itemindent
              \parsep \z@
            \fi
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand{\theenumiv}{\arabic{enumiv}}}%
      \if@openbib
        \renewcommand{\newblock}{\par}
      \else
        \renewcommand{\newblock}{\hskip .11em \@plus.33em \@minus.07em}%
      \fi
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand{\newblock}{}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INDEX
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\@makessectionhead{\indexname}]%
                \@mkboth{\uppercase{\indexname}}%
                        {\uppercase{\indexname}}%
                \thispagestyle{class}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand{\@idxitem}  {\par\hangindent 40\p@}
\newcommand{\subitem}   {\par\hangindent 40\p@ \hspace*{20\p@}}
\newcommand{\subsubitem}{\par\hangindent 40\p@ \hspace*{30\p@}}
\newcommand{\indexspace}{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand{\footnoterule}{%
  \kern-3\p@
  \hrule width .4\columnwidth
  \kern 2.6\p@}
\long\def\@makefntext#1{%
    \parindent 1em%
    \noindent
    \hbox to 1.8em{\hss\@makefnmark}#1}
\newcommand{\contentsname}{Contents}
\newcommand{\listfigurename}{List of Figures}
\newcommand{\listtablename}{List of Tables}
\newcommand{\refname}{References}
\newcommand{\bibname}{Bibliography}
\newcommand{\indexname}{Index}
\newcommand{\figurename}{Figure}
\newcommand{\tablename}{Table}
\newcommand{\partname}{Part}
\newcommand{\sectionname}{Section}
\newcommand{\annexname}{Annex}
\newcommand{\appendixname}{Appendix}
\newcommand{\executivename}{Executive Summary}
\newcommand{\abstractname}{Abstract}
\newcommand{\keywordsname}{Keywords}
\newcommand{\today}{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{class}
\pagenumbering{arabic}      % Arabic page numbers
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi
\endinput
%% 
%% End of file `saclantc.cls'.
