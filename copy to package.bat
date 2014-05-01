md package
copy package_readme. package\README
robocopy . package /xf package_readme README .git* *.bat *.pyc *.zip
robocopy "sphinx/_build/html" "package/doc/html" /s /xf .buildinfo objects.inv
robocopy "example" "package/example" /s /xf *.rst *.py
