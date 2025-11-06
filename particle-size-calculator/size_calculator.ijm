//GGE particle size macro for ImageJ to replace GGE7A3 on NIHImage.
imageID = getImageID();

main();

//Global vars, will transition away from global vars once working properly.

var stdWeights;
var stdRfVals;
var imageWidth;
var imageHeight;
var xnum;
var c = 0;
var z = 0;
var t = 0;
var cubicCoeffArray;
var quarticCoeffArray;
var imageID;
var originXVal;
var laneLength;
var xValsCurrLaneIndex;
var yValsCurrLane;
var bins;



// Main macro flow
function main() {
    //print('main');
    initialize();
    getStandards();
    quarticCoeffArray = calcStdCurveQuartic(stdRfVals, stdWeights);
    //cubicRegression(stdRfVals, stdWeights);
    print('Regression coefficients:');
    Array.print(quarticCoeffArray);
    getLanes();
}

// Prompts user to crop gel
function initialize() {
    //print('initialize');
    if (nImages == 0) {
        print('Error: Open an image before running Macro.');
        selectWindow('Log');
        exit();
    }
    if (nImages > 1) {
        print('Error: Multiple images detected. Have one image open before running Macro.');
        selectWindow('Log');
        exit();
    }
    imageWidth = 0;
    imageHeight = 0;
    xnum = getTitle();

    run('Set Measurements...', 'invert redirect=None decimal=3');
    run('Gel Analyzer Options...', 'vertical=1 horizontal=1 label');
    run('Clear Results');
    run("Overlay Options...", "stroke=red width=1 apply");
    print('');
    print('Analysis of ' + xnum);

    setTool('Rectangle');

    getDimensions(imageWidth, imageHeight, c, z, t);
    makeRectangle(0.3*imageWidth, 0.18*imageHeight, 0.36*imageWidth, 0.53*imageHeight);
    waitForUser('Adjust rectangle to bound gel. \nPressing OK will crop gel and rotate 90 degrees left.');

    if (selectionType() != 0) {
        print('Error: Selection must be rectangular.');
        selectWindow('Log');
        exit();
    }

    run('Crop');
    wait(100);
    run('Rotate 90 Degrees Left');

}

function getStandards() {
    //print('getStandards');
    stdWeights = newArray(0);
    stdRfVals = newArray(0);
    moreStandards = true;
    while (moreStandards) {
        stdWeightsTempAndStdRfValsTemp = setStandard();
        arrayLength = stdWeightsTempAndStdRfValsTemp.length;
        stdWeightsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, 0, (arrayLength/2));
        stdRfValsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, (arrayLength/2), arrayLength);
        stdWeights = Array.concat(stdWeights, stdWeightsToAdd);
        stdRfVals = Array.concat(stdRfVals, stdRfValsToAdd);

        moreStandards = false;
        wait(300);
        moreStandards = getBoolean('Add more standard lanes?');
    }
    Array.print(stdWeights);
    Array.print(stdRfVals);
}

function setStandard() {
    //print('setStandards');
    selectImage(imageID);
    numberOfStandards = getNumber('enter number of standards', 4);
    stdWeightsTemp = newArray(numberOfStandards);
    for(i = 0; i < numberOfStandards; i++){
        stdIndexDisplay = i + 1;
        stdWeightsTemp[i] = getNumber('enter weight of standard ' + stdIndexDisplay, 0);
    }

    setTool('Rectangle');

    selectImage(imageID);
    getDimensions(standardLaneWidth, standardLaneHeight, c, z, t);

    makeRectangle(0, (standardLaneHeight/2)-10 , standardLaneWidth, 20);
    waitForUser('Adjust rectangle to span lane with standards. \nIt can be quite thin as long as it contains some part of the lane.');
    
    run('Select First Lane');
    stdRfValsTemp = getRfValsFromLaneLineCursor('standard');
    //Array.print(stdRfValsTemp);
    stdWeightsTempAndStdRfValsTemp = Array.concat(stdWeightsTemp, stdRfValsTemp);
    //Array.print(stdWeightsTempAndStdRfValsTemp);
    if (2 * numberOfStandards != stdWeightsTempAndStdRfValsTemp.length) {
        print('Std weights or Rf vals not recorded properly.');
    }
    return stdWeightsTempAndStdRfValsTemp;
}


function getRfValsFromLaneLineCursor(peakType) {
    if (selectionType() != 0) {
        print('Error: Selection must be rectangular.');
        selectWindow('Log');
        exit();
    }
    waitForUser('When lane plot is displayed:\n1. Select the origin peak\n2. Select ' + peakType + ' peaks\n3. Press space bar when finished');
    run('Plot Profile');
    run('Remove Overlay');
    setTool('multi-point');
    run('Clear Results');
    
    Plot.getValues(xValsCurrLane, yValsCurrLane);
    xValsCurrLaneIndex = divideArrayByStep(xValsCurrLane, xValsCurrLane[1]-xValsCurrLane[0]);

    getDimensions(w, h, c, z, f);
    
    while(true) {

        getCursorLoc(x, y, z, m);
        Overlay.drawLine(x, 0, x, h);
        Overlay.drawLine(0, y, w, y);
        Overlay.show();
        wait(1);
        Overlay.remove();
        if (isKeyDown('space')) {
            run('Measure');
            wait(300);
            xValsAndOrigin = divideArrayByStep(getAllResults('X'), xValsCurrLane[1]-xValsCurrLane[0]);

            originXVal = xValsAndOrigin[0];
            laneLength = xValsCurrLaneIndex.length-originXVal;
            xValsCurrLaneIndex = Array.slice(xValsCurrLaneIndex, originXVal, xValsCurrLaneIndex.length);

            for (i = 0; i < xValsCurrLaneIndex.length; i ++) {
                xValsCurrLaneIndex[i] = xValsCurrLaneIndex[i] - originXVal;
            }

            yValsCurrLane = Array.slice(yValsCurrLane, originXVal, yValsCurrLane.length);

            xVals = Array.slice(xValsAndOrigin, 1);

            for (i = 0; i < xVals.length; i ++) {
                xVals[i] = xVals[i] - originXVal;
            }
            //print('xVals');
            //Array.print(xVals);
            // std xvals are ok here
            RfVals = calcRfVals(xVals, laneLength);
            return RfVals;
        }                
    }
    print('past while loop');
    
}

function getLanes() {
    //print('getLanes');
    moreLanes = true;
    while (moreLanes) {
        quantifyLane();
        moreLanes = getBoolean('Analyze more lanes?');
    }
}

function quantifyLane() {
    //print('quantifyLane');
    if (stdRfVals.length == 0 || stdWeights.length == 0){
        exit('Standards not set.');
    }

    selectImage(imageID);

    waitForUser('Drag rectangle to sample lane and press OK.');

    run('Select First Lane');

    laneRfVals = getRfValsFromLaneLineCursor('unknown');

    //laneMolecularWeightsCalc = newArray(laneRfVals.length);
    // for (i = 0; i < laneRfVals.length; i++) {
    //     x = laneRfVals[i];
    //     laneMolecularWeightsCalc[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*x + cubicCoeffArray[2]*x*x + cubicCoeffArray[3]*x*x*x;
    // }
    laneMolecularWeightsCalc = predictQuartic(laneRfVals, quarticCoeffArray);
    displayValueArray = inverseLog10ArrayDisabled(laneMolecularWeightsCalc);
    print('Calculated Diameters:');
    Array.print(displayValueArray);
    print('');
    print('LDL Bins:');
    quantBins();

}

function quantBins() {
    //print('quantBins');
    binPxValues = calcPxFromBins();
    baselineY = getBackgroundConc();
    Array.getStatistics(yValsCurrLane, min, max, mean, stdDev);
    //print(min, max, mean, stdDev);
    //Array.print(yValsCurrLane);
    for (i = 0; i < yValsCurrLane.length; i++) {
    yValsCurrLane[i] = yValsCurrLane[i] - baselineY;
    }
    binSums = newArray(binPxValues.length-1);
    for (i = 0; i < binPxValues.length-1; i ++) {
        binSums[i] = sumSingleBin(yValsCurrLane, binPxValues[i], binPxValues[i+1]);
    }
    Array.getStatistics(binSums, min, max, mean, stdDev);
    binSumsTotal = mean*binSums.length;
    for (i = 0; i < binSums.length; i++) {
        binSums[i] = (binSums[i]/binSumsTotal)*100;
        print(bins[i] + '-' + bins[i+1] + ': ' + binSums[i] + '%');
    }
    
    //Array.print(binSums);

}


function calcPxFromBins() {
    //print('calcPxFromBins');
    bins = newArray(
    375, 339, 321, 315, 309, 303, 297, 291, 285, 272, 265, 256, 247, 242, 233, 220
    );
    binCount = bins.length; 
    everyRfValue = calcRfVals(xValsCurrLaneIndex, laneLength);
    // everyLogMW = newArray(xValsCurrLaneIndex.length);
    // for (i = 0; i < laneLength; i ++) {
    //         everyLogMW[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*everyRfValue[i] + cubicCoeffArray[2]*everyRfValue[i]*everyRfValue[i] + cubicCoeffArray[3]*everyRfValue[i]*everyRfValue[i]*everyRfValue[i];
    // }
    everyLogMW = predictQuartic(everyRfValue, quarticCoeffArray);
    everyMW = inverseLog10ArrayDisabled(everyLogMW);
    binPxValues = newArray(binCount);

    for (j = 0; j < bins.length; j ++) {
        everyMWCopy = Array.copy(everyMW);
        for (i = 0; i < laneLength; i ++) {
            everyMWCopy[i] = abs(everyMWCopy[i] - bins[j]);   
        }
        Array.getStatistics(everyMWCopy, min, max, mean, stdDev);

        for (h = 0; h < laneLength; h ++) {
            if (abs(everyMWCopy[h] - min)<1e-6) {
                binPxValues[j] = h;
            }
        }
    }
    //Array.print(binPxValues);
    return binPxValues;
}

function sumSingleBin(yValArray, binLowerBound, binUpperBound) {
    //print('sumSingleBin');
    slicedArray = Array.slice(yValArray, binLowerBound, binUpperBound);
    //Array.print(slicedArray);
    Array.getStatistics(slicedArray, min, max, mean, stdDev);
    //print(min, max, mean, stdDev);
    //Array.print(slicedArray);
    return slicedArray.length * mean;
}

function cubicRegression(stdRfValues, stdWeights) {
    //print('cubicRegression');
    log10StdWeights = log10ArrayDisabled(stdWeights);
    cubicCoeffArray = calcStdCurveCubic(stdRfVals, log10StdWeights);
    
    //Array.print(cubicCoeffArray);
    return cubicCoeffArray;
}

function getBackgroundConc() {
    //print('getBackgroundConc');

    waitForUser('Select point that reflects baseline y-value of the LDL range. Press space bar when finished. \nPress OK before returning to lane plot.');
    run('Remove Overlay');
    setTool('multi-point');
    
    getDimensions(w, h, c, z, f);
    
    while(true) {

        getCursorLoc(x, y, z, m);
        Overlay.drawLine(x, 0, x, h);
        Overlay.drawLine(0, y, w, y);
        Overlay.show();
        
        if (isKeyDown('space')) {
            resultsCountBefore = nResults;
            run('Clear Results');
            run('Measure');
            if (selectionType() != 10 || nResults - resultsCountBefore != 1) {
                print('Error: Invalid selection. Selection must be a single point.');
                selectWindow('Log');
                exit();
            }
            wait(300);
            return getResult('Y');
        }

        wait(1);
        Overlay.remove();
                    
    }
}



// -- Utils -- //

// Rounds number to given decimal step
function roundToStep(number, step) {
    return round(number/step)*step;
}

// roundToStep for array
function roundArrayToStep(array, step) {
    roundedArray = newArray(array.length);
    for (i = 0; i < array.length; i ++) {
        roundedArray[i] = roundToStep(array[i], step);
    }
    return roundedArray;
}

// Returns whole multiples of step in number
function divideByStep(number, step) {
    return round(number/step);
}

// divideByStep for array
function divideArrayByStep(array, step) {
    dividedArray = newArray(array.length);
    for (i = 0; i < array.length; i ++) {
        dividedArray[i] = divideByStep(array[i], step);
    }
    return dividedArray;
}

// Retrieves all results currently in Results window
function getAllResults(column) {
    array = newArray(nResults);
    results = newArray(array.length);
    for (i = 0; i < results.length; i ++) {
        results[i] = getResult(column, i);
    }
    return results;
}

// Calculates Rf vaules given an array of x values and total lane length
function calcRfVals(xVals, laneLength) {
    //print('calcRfVals');
    result = newArray(xVals.length);
    for (i = 0; i < xVals.length; i ++) {
        result[i] = xVals[i]/laneLength;
    }
    //print('Rf');
    //Array.print(result);
    return result;
}

// Calculates log base 10 for every value of an array
function log10ArrayDisabled(array) {
    //print('log10Array');
    // logValues = newArray(array.length);
    // for (i = 0; i < array.length; i++){
    //     logValues[i] = log(array[i])/log(10);
    // }
    return array;
}

// Calculates inverse log 10 for every value of an array
function inverseLog10ArrayDisabled(array) {
    //print('inverseLog10Array');
    // inverseLogValues = newArray(array.length);
    // for (i = 0; i < array.length; i++){
    //     inverseLogValues[i] = pow(10, array[i]);
    // }
    return array;
}

// Performs OLS regression for arrays of x values and y values
function OLSRegression(xValsForRegression, yValsForRegression) {
    //print('OLSRegression');
    xLength = xValsForRegression.length;
    yLength = yValsForRegression.length;

    xTotal = 0;
    yTotal = 0;

    for(i = 0; i < xLength; i++) {
        xTotal += xValsForRegression[i];
        yTotal += yValsForRegression[i];
    }
    xMean = xTotal/xLength;
    yMean = yTotal/yLength;

    numerator = 0;
    denominator = 0;

    for(i = 0; i < xLength; i++) {
        dx = xValsForRegression[i] - xMean;
        dy = yValsForRegression[i] - yMean;
        numerator += dx * dy;
        denominator += dx * dx;
    }

    regressionSlope = numerator/denominator;

    regressionIntercept = yMean - regressionSlope * xMean;

    return newArray(regressionSlope, regressionIntercept);
}

// Calculates cubic regression coefficients given x values and y values
function calcStdCurveCubic(xVals, yVals) {
    //print('calcStdCurveCubic');
    n = xVals.length;

    sumX = newArray(7);
    for (i = 0; i <= 6; i++) sumX[i] = 0;
    sumXY = newArray(4);
    for (i = 0; i <= 3; i++) sumXY[i] = 0;

    for (i = 0; i < n; i++) {
        x = xVals[i];
        y = yVals[i];
        powX = newArray(7);
        powX[0] = 1;
        for (j = 1; j <= 6; j++) powX[j] = powX[j - 1] * x;

        for (j = 0; j <= 6; j++) sumX[j] += powX[j];
        for (j = 0; j <= 3; j++) sumXY[j] += powX[j] * y;
    }

    // Flattened 4x4 matrix A
    A = newArray(16);
    for (i = 0; i <= 3; i++) {
        for (j = 0; j <= 3; j++) {
            A[i * 4 + j] = sumX[i + j];
        }
    }

    B = sumXY;

    coeffs = gaussJordanFlat4x4(A, B);
    return coeffs; // [a, b, c, d]
}

// Gauss Jordan reduction of 4x4 matrix
function gaussJordanFlat4x4(A, B) {
    //print('gaussJordanFlat4x4');
    n = 4;

    for (i = 0; i < n; i++) {
        // Find non-zero pivot
        if (A[i * 4 + i] == 0) {
            for (j = i + 1; j < n; j++) {
                if (A[j * 4 + i] != 0) {
                    // Swap rows in A
                    for (k = 0; k < n; k++) {
                        temp = A[i * 4 + k];
                        A[i * 4 + k] = A[j * 4 + k];
                        A[j * 4 + k] = temp;
                    }
                    // Swap B
                    tmpB = B[i];
                    B[i] = B[j];
                    B[j] = tmpB;
                    break;
                }
            }
        }

        // Normalize row
        factor = A[i * 4 + i];
        for (j = 0; j < n; j++) A[i * 4 + j] /= factor;
        B[i] /= factor;

        // Eliminate other rows
        for (j = 0; j < n; j++) {
            if (j != i) {
                factor = A[j * 4 + i];
                for (k = 0; k < n; k++) {
                    A[j * 4 + k] -= factor * A[i * 4 + k];
                }
                B[j] -= factor * B[i];
            }
        }
    }

    return B;
}

// Calculates quartic regression coefficients 
function calcStdCurveQuartic(xVals, yVals) {
    //print('calcStdCurveQuartic');
    n = xVals.length;

    sumX = newArray(9); // sums of x^0 ... x^8
    for (i = 0; i <= 8; i++) sumX[i] = 0;
    sumXY = newArray(5); // sums of x^0*y ... x^4*y
    for (i = 0; i <= 4; i++) sumXY[i] = 0;

    for (i = 0; i < n; i++) {
        x = xVals[i];
        y = yVals[i];
        powX = newArray(9);
        powX[0] = 1;
        for (j = 1; j <= 8; j++) powX[j] = powX[j - 1] * x;

        for (j = 0; j <= 8; j++) sumX[j] += powX[j];
        for (j = 0; j <= 4; j++) sumXY[j] += powX[j] * y;
    }

    // Flattened 5x5 matrix A
    A = newArray(25);
    for (i = 0; i <= 4; i++) {
        for (j = 0; j <= 4; j++) {
            A[i * 5 + j] = sumX[i + j];
        }
    }

    B = sumXY;

    coeffs = gaussJordanFlat5x5(A, B);
    return coeffs; // [a, b, c, d, e]
}

// Gauss Jordan reduction of 5x5 matrix
function gaussJordanFlat5x5(A, B) {
    //print('gaussJordanFlat5x5');
    n = 5;

    for (i = 0; i < n; i++) {
        // Find non-zero pivot
        if (A[i * n + i] == 0) {
            for (j = i + 1; j < n; j++) {
                if (A[j * n + i] != 0) {
                    // Swap rows in A
                    for (k = 0; k < n; k++) {
                        temp = A[i * n + k];
                        A[i * n + k] = A[j * n + k];
                        A[j * n + k] = temp;
                    }
                    // Swap B
                    tmpB = B[i];
                    B[i] = B[j];
                    B[j] = tmpB;
                    break;
                }
            }
        }

        // Normalize row
        factor = A[i * n + i];
        for (j = 0; j < n; j++) A[i * n + j] /= factor;
        B[i] /= factor;

        // Eliminate other rows
        for (j = 0; j < n; j++) {
            if (j != i) {
                factor = A[j * n + i];
                for (k = 0; k < n; k++) {
                    A[j * n + k] -= factor * A[i * n + k];
                }
                B[j] -= factor * B[i];
            }
        }
    }

    return B;
}

// Estimates y value from x value based on quartic regression
function predictQuartic(xVals, quarticCoeffArray) {
    result = newArray(xVals.length);
    for (i = 0; i < xVals.length; i ++) {
        result[i] = quarticCoeffArray[0] + 
                    quarticCoeffArray[1]*xVals[i] + 
                    quarticCoeffArray[2]*xVals[i]*xVals[i] + 
                    quarticCoeffArray[3]*xVals[i]*xVals[i]*xVals[i] +
                    quarticCoeffArray[4]*xVals[i]*xVals[i]*xVals[i]*xVals[i];
    }
    return result;
}

// Connects cursor to horizontal line
function addHorizontalLineToCursor() {
    getDimensions(w, h, c, z, f);
    while(true) {
        getCursorLoc(x, y, z, m);
        Overlay.drawLine(0, y, w, y);
        Overlay.show();
        wait(5);
        Overlay.remove();
}
}


// Connects cursor to vertical line
function addVerticalLineToCursor() {
    getDimensions(w, h, c, z, f);
    while(true) {
        getCursorLoc(x, y, z, m);
        Overlay.drawLine(x, 0, x, h);
        Overlay.show();
        wait(5);
        Overlay.remove();
    }
}