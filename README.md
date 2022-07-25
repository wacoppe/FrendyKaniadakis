# FrendyKaniadakis
:radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive: :radioactive:  

## Introduction

Hello, 


Welcome to our GitHub repository. We are members of the reactors physics laboratory at the Federal University of Rio de Janeiro ([UFRJ/COPPE](https://www.coppe.ufrj.br/en)) and this is our modifications in the nuclear data processing code called FRENDY, developed by the Japan Atomic Energy Agency ([JAEA](https://www.jaea.go.jp/english/)). Our primary goal is to generate deformed neutron cross-sections using the Kaniadakis statistics.

(_This topic will be continually improved with the involved theory_)

## How to run

In order to TEST these modifications, it is **mandatory** to install the FRENDY code in your machine. The following website explains how to do that:
https://rpg.jaea.go.jp/main/en/program_frendy/

These are the files that will be modified inside FRENDY
in order to calculating deformed cross-sections using the Kaniadakis distribution

#### Utilized folders inside FRENDY:

* _frendy_XXXXXXXX\tests\ReconResonance_ (File "_ResonanceXSCalculatorTest.cpp_" and Folder "_for_test_" )

* _frendy_XXXXXXXX\frendy\ReconResonance_ (Files "_ResonanceXSCalculator.cpp_" and "_ResonanceXSCalculator.hpp_")


#### Folder to see the results files:

* _frendy_XXXXXXXX\tests\ReconResonance\comp_njoy_

#### Folder to start the compilation (Windows, via 'wsl' command):

* _frendy_XXXXXXXX\tests\ReconResonance_


## The Fadeeva Package

In order to calculate some critical Faddeeva functions inside our analytical solution for the deformed Doppler broadening function, we implemented the [Faddeeva Package](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package), developed by the MIT. Consequently, the modifications will only work with the Faddeeva files in both aforementioned folders. We thank Professor [Steven G. Johnson](https://math.mit.edu/~stevenj/) for his incredible work.


## About us

Programa de Engenharia Nuclear ([PEN - COPPE/UFRJ](http://www.con.ufrj.br/)) (Brazil)    
Editors: Willian Abreu (wabreu@coppe.ufrj.br) / João Maciel (joaomarcio.maciel@poli.ufrj.br)    
Supervisor: Prof. Aquilino Senra Martinez  

## Acknowledgments

Willian Abreu and his supervisor, Professor Aquilino Senra Martinez, thank [**FAPERJ**](https://www.faperj.br/) for the conceded grant (Pós-doutorado Nota 10). Professor Aquilino Senra Martinez also thank the National Council for Scientific and Technological Development ([CNPq](https://www.gov.br/cnpq/pt-br)).

| Support/funding:     |    |       |   |
| :---         |     :---:      |          ---: | :---    |
| <img src="http://www.con.ufrj.br/wp-content/uploads/2015/07/logo.gif" width="200">   | <img src="https://www.faperj.br/downloads/logomarcas/logo.jpg" width="200" style="text-align:center">     | <img src="https://www.gov.br/cnpq/pt-br/canais_atendimento/identidade-visual/CNPq_v2017_rgb.jpg" width="200" > | <img src="https://upload.wikimedia.org/wikipedia/pt/1/1e/Logo_COPPE_-_UFRJ.jpg" width="200" style="text-align:center">




