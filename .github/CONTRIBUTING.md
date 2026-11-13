# Cómo contribuir / Contributing

**[Español](#español) | [English](#english)**

## Español

¡Gracias por tu interés en GeneExprLupusMx! Puedes escribir en español o en inglés.

### Formas de contribuir

- **Traducciones a otros idiomas.** Las lecciones están en español e inglés. Si quieres traducirlas a otro idioma (por ejemplo, portugués), abre un *issue* para organizarnos antes de empezar.
- **Revisión de los materiales.** Si encuentras un error, una explicación poco clara o un enlace roto en las lecciones, el README o la ayuda de las funciones, avísame en un *issue*.
- **Sugerencias de gráficas y análisis.** ¿Se te ocurre otra gráfica u otro análisis que se pueda hacer con estos datos? Cuéntame tu idea en un *issue*.
- **Comparte tus ideas o tu código.** Si ya hiciste un análisis con el paquete, puedes compartirlo en un *issue* o en un *pull request*.
- **Reportar errores.** Si una función no hace lo que esperas, abre un *issue* con un ejemplo mínimo reproducible ([reprex](https://reprex.tidyverse.org/)).

### Dónde está cada cosa

- **Lecciones:** `vignettes/articles/`, un archivo por idioma (`*_es.Rmd` y `*_en.Rmd`). Si cambias una lección, cambia también su versión en el otro idioma (o avísame y yo la actualizo).
- **Funciones y documentación de los datos:** `R/`. La ayuda (`man/`) se genera con [roxygen2](https://roxygen2.r-lib.org/): edita los comentarios `#'` en `R/`, no los archivos `.Rd`.
- **Pruebas:** `tests/testthat/`.

### *Pull requests*

1. Haz un *fork* y clónalo, por ejemplo con `usethis::create_from_github("EveliaCoss/GeneExprLupusMx", fork = TRUE)`.
2. Instala todas las dependencias, incluidas las de `Suggests`, con `pak::local_install_dev_deps()`, y crea una rama con `usethis::pr_init("descripcion-breve")`.
3. Haz tus cambios. Si cambias una función, agrega o actualiza sus pruebas y comprueba que `devtools::check()` pase sin errores.
4. Si el cambio es visible para quien usa el paquete, agrega una línea al inicio de `NEWS.md`.
5. Envía el *pull request* con `usethis::pr_push()`.

Para cambios grandes, abre primero un *issue* para que lo platiquemos.

### Código de conducta

Este paquete tiene un [Código de Conducta](CODE_OF_CONDUCT.md), basado en el Contributor Covenant. Al participar, aceptas cumplirlo.

## English

Thank you for your interest in GeneExprLupusMx! You can write in Spanish or in English.

### Ways to contribute

- **Translations into other languages.** The lessons are in Spanish and English. If you would like to translate them into another language (for example, Portuguese), open an issue so we can organize before you start.
- **Review of the materials.** If you find a mistake, an unclear explanation or a broken link in the lessons, the README or the function help, let me know in an issue.
- **Suggestions of plots and analyses.** Can you think of another plot or analysis that could be done with these data? Tell me about your idea in an issue.
- **Share your ideas or your code.** If you already did an analysis with the package, you can share it in an issue or in a pull request.
- **Report bugs.** If a function does not do what you expect, open an issue with a minimal reproducible example ([reprex](https://reprex.tidyverse.org/)).

### Where things are

- **Lessons:** `vignettes/articles/`, one file per language (`*_es.Rmd` and `*_en.Rmd`). If you change a lesson, also change its version in the other language (or let me know and I will update it).
- **Functions and data documentation:** `R/`. The help files (`man/`) are generated with [roxygen2](https://roxygen2.r-lib.org/): edit the `#'` comments in `R/`, not the `.Rd` files.
- **Tests:** `tests/testthat/`.

### Pull requests

1. Fork and clone the repository, for example with `usethis::create_from_github("EveliaCoss/GeneExprLupusMx", fork = TRUE)`.
2. Install all the dependencies, including those in `Suggests`, with `pak::local_install_dev_deps()`, and create a branch with `usethis::pr_init("brief-description")`.
3. Make your changes. If you change a function, add or update its tests and check that `devtools::check()` passes without errors.
4. If the change is visible to package users, add a bullet at the top of `NEWS.md`.
5. Submit the pull request with `usethis::pr_push()`.

For big changes, please open an issue first so we can discuss it.

### Code of Conduct

This package has a [Code of Conduct](CODE_OF_CONDUCT.md), based on the Contributor Covenant. By participating, you agree to abide by its terms.
