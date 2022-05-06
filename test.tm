<TeXmacs|2.1.1>

<style|<tuple|generic|maxima>>

<\body>
  <\session|maxima|default>
    <\output>
      \;

      Maxima 5.46.0 https://maxima.sourceforge.io

      using Lisp SBCL 2.1.9

      Distributed under the GNU Public License. See the file COPYING.

      Dedicated to the memory of William Schelter.

      The function bug_report() provides bug reporting information.
    </output>

    <\input>
      <with|color|red|(<with|math-font-family|rm|%i>1) >
    <|input>
      B : [24.0806, 11.3630, 31.1384, 66.1825, 26.2870]$
    </input>

    <\input>
      <with|color|red|(<with|math-font-family|rm|%i>2) >
    <|input>
      A : matrix(

      [48.6760, 1.62932, 1.88423, 1.63136, 1.58985],

      [1.92203, 56.6858, 1.21301, 1.37236, 1.06716],

      [1.66020, 1.38067, 12.2637, 1.41381, 1.95707],

      [1.90243, 1.52317, 1.16776, 46.1008, 1.56883],

      [1.64947, 1.50724, 1.69159, 1.60908, 12.1297]

      )$
    </input>

    <\unfolded-io>
      <with|color|red|(<with|math-font-family|rm|%i>3) >
    <|unfolded-io>
      float(invert(A).B)
    <|unfolded-io>
      <math|<with|math-display|true|<text|<with|font-family|tt|color|red|(<with|math-font-family|rm|%o3>)
      >><matrix|<tformat|<table|<row|<cell|0.31383834993401>>|<row|<cell|0.08263252251914907>>|<row|<cell|2.072610817993732>>|<row|<cell|1.311230448655671>>|<row|<cell|1.651228273556617>>>>>>>
    </unfolded-io>

    <\input>
      <with|color|red|(<with|math-font-family|rm|%i>4) >
    <|input>
      \;
    </input>
  </session>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>