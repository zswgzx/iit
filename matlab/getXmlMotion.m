function [interlaceAngle,interlaceTrans,gradientAngle,gradientTrans] = getXmlGradientMotion(filename)

% see http://blogs.mathworks.com/community/2010/11/01/xml-and-matlab-navigating-a-tree/

gradMotion=xmlread(filename);

% get the xpath mechanism into the workspace
import javax.xml.xpath.*
factory = XPathFactory.newInstance;
xpath = factory.newXPath;

numVol=84;
interlaceAngle=zeros(numVol,3);
interlaceTrans=zeros(numVol,3);
gradientAngle=zeros(numVol,3);
gradientTrans=zeros(numVol,3);

for i=0:numVol-1
    gradient=sprintf('\''gradient_%04d\''',i);
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceAngleX'']/green'];
    expression = xpath.compile(expr);
    interlaceAngle(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(interlaceAngle(i+1,1))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceAngleX'']/red'];
        expression = xpath.compile(expr);
        interlaceAngle(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceAngleY'']/green'];
    expression = xpath.compile(expr);
    interlaceAngle(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(interlaceAngle(i+1,2))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceAngleY'']/red'];
        expression = xpath.compile(expr);
        interlaceAngle(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceAngleZ'']/green'];
    expression = xpath.compile(expr);
    interlaceAngle(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(interlaceAngle(i+1,3))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceAngleZ'']/red'];
        expression = xpath.compile(expr);
        interlaceAngle(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceTranslationX'']/green'];
    expression = xpath.compile(expr);
    interlaceTrans(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(interlaceTrans(i+1,1))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceTranslationX'']/red'];
        expression = xpath.compile(expr);
        interlaceTrans(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceTranslationY'']/green'];
    expression = xpath.compile(expr);
    interlaceTrans(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(interlaceTrans(i+1,2))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceTranslationY'']/red'];
        expression = xpath.compile(expr);
        interlaceTrans(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceTranslationZ'']/green'];
    expression = xpath.compile(expr);
    interlaceTrans(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(interlaceTrans(i+1,3))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''InterlaceWiseCheck'']/entry[@parameter=''InterlaceTranslationZ'']/red'];
        expression = xpath.compile(expr);
        interlaceTrans(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end    
    
    
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientAngleX'']/green'];
    expression = xpath.compile(expr);
    gradientAngle(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(gradientAngle(i+1,1))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientAngleX'']/red'];
        expression = xpath.compile(expr);
        gradientAngle(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientAngleY'']/green'];
    expression = xpath.compile(expr);
    gradientAngle(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(gradientAngle(i+1,2))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientAngleY'']/red'];
        expression = xpath.compile(expr);
        gradientAngle(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientAngleZ'']/green'];
    expression = xpath.compile(expr);
    gradientAngle(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(gradientAngle(i+1,3))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientAngleZ'']/red'];
        expression = xpath.compile(expr);
        gradientAngle(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientTranslationX'']/green'];
    expression = xpath.compile(expr);
    gradientTrans(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(gradientTrans(i+1,1))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientTranslationX'']/red'];
        expression = xpath.compile(expr);
        gradientTrans(i+1,1)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientTranslationY'']/green'];
    expression = xpath.compile(expr);
    gradientTrans(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(gradientTrans(i+1,2))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientTranslationY'']/red'];
        expression = xpath.compile(expr);
        gradientTrans(i+1,2)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end
    
    expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientTranslationZ'']/green'];
    expression = xpath.compile(expr);
    gradientTrans(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    if isnan(gradientTrans(i+1,3))==1
        expr=['/QCResultSettings/entry[@parameter=''DWI Check'']/entry[@parameter=',gradient,']/entry[@parameter=''GradientWiseCheck'']/entry[@parameter=''GradientTranslationZ'']/red'];
        expression = xpath.compile(expr);
        gradientTrans(i+1,3)=expression.evaluate(gradMotion,XPathConstants.NUMBER);
    end    
end
