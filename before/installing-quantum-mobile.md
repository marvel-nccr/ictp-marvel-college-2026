# How to install Quantum Mobile

## Installation instructions

Get Quantum Mobile running on your computer in three simple steps:

 1. Download the appropriate image for your architecture:
    
    | Architecture | Link | 
    |---|---|
    | Intel x86/AMD64     | [Google Drive](https://drive.google.com/file/d/12hFWtAFn83wqAwF9H6ClXbb9adOQfpeb/view?usp=sharing) |
    | Apple Silicon/ARM64 | [Google Drive](https://drive.google.com/file/d/1Xqogj1__GRIq2heB92s_2ff3QUF1nAZG/view?usp=sharing) |

 2. Install Virtual Box 7.2.6 or later (see <https://www.virtualbox.org>)
 3. Import the ``.ova`` file into Virtualbox: ``File`` ➔ ``Import Appliance``

## Troubleshooting
Solutions to common problems can be found [here](https://quantum-mobile.readthedocs.io/en/latest/users/troubleshoot.html). If that page does not answer your questions, feel free to ask questions in ``#quantum-mobile`` on the ICTP-MARVEL slack

## Tips
By default, Quantum Mobile uses a US keyboard layout. To change it...
- open the ubuntu settings within the VM (click on the ubuntu logo in the lower left corner, search for ``Settings``, and then click on it to launch)
- go to ``Keyboard`` ➔ ``Input sources`` ➔ ``Add Input Source``
- click on the ``⋮`` and add the desired language
- Once added there you should be able to select the desired keyboard profile in a drop-down menu on the upper right corner of the screen.