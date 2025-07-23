window.MathJax = {
  loader: {load: ['[tex]/upgreek']},
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    packages: {'[+]': ['upgreek']}
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};
  
document$.subscribe(() => {
  //MathJax.typesetPromise()
})