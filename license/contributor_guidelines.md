Document detailing what is expected of contributors, to be sent to new
contributors when they get commit access.

**NOTE**: this is only a draft, to be finalized after further discussion.
**NOTE**: PDF files mentioned below can be generated using `make`

# Social aspects

## Code of conduct

Your interactions with the community should follow the code described in `license/CODE_OF_CONDUCT.md`

## Contactability

Add your contact details to `license/developers.txt`.

When you make a commit using Git, your email address is stored with the commit metadata.
Make sure the email address you use is one you regularly monitor.
In case there are any problems with your commit, you may be contacted by other developers.
To change the email address used for your commits, see
<https://stackoverflow.com/questions/37805621/change-email-address-in-git/37805844#37805844>.

## Mailing lists

Major announcements regarding the code will be made on
<https://groups.google.com/g/pencil-code-discuss>;
it is recommended that you subscribe to this mailing list.
Changes to the Python module can be discussed on
<https://groups.google.com/g/pencil-code-python>

# Recommended development style

## Version control
**TODO**: under discussion

To effectively use Git, see `doc/git/git-best-practises.pdf`

## Coding style

### Fortran

See section 9.1 ('The PENCIL CODE coding standard') of the manual (`doc/manual.pdf`)

### Python
- Four spaces (not tabs) for indentation

## What to do before committing changes

### Do I need to consult anyone?

If you change the reference data for any of the samples (the `reference.out`
files in the corresponding directories), you *must* contact the maintainer of
that particular sample (listed in the README file) and get their permission.
In your commit, mention which maintainer approved your change.

Other than the above condition, changes are welcome as long as they do not break
existing functionality (see the next section, 'Tests').
If you think your changes may disrupt other users' usage of the code, we
recommend soliciting feedback during the Pencil Office Hours or on the mailing
list *before* you add your changes to the master branch.

### Tests

See section 10 ('Testing the code') of the manual.
If you make changes to the postprocessing (e.g. Python) modules, make sure you
also run the associated tests (`pc_auto-test` mainly tests the core Fortran
part of the code).
