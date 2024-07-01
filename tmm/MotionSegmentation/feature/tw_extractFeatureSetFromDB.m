function [ featureSet, db ] = tw_extractFeatureSetFromDB( db, featureType, frame_offsets)
%TW_EXTRACTFEATURESETFROMDB Summary of this function goes here
%   Detailed explanation goes here

    switch nargin
        case 2
            frames = 1:db.nrOfFrames;
        case 3
            frames = 1:db.nrOfFrames;
        otherwise
            error('\nInvalid number of arguments!\n');
    end

    switch lower(featureType)
        case 'e15' 
            featureSet = tw_extractFeatureSetFromDB_local(db,frames,'pos',[4,9,17,21,28]);
        %----------------------------------------------------------------------
        case 'e27' 
            featureSet = tw_extractFeatureSetFromDB_local(db,frames,'pos',[3,4,8,9,17,19,21,26,28]);
        %----------------------------------------------------------------------
        case 'e30' 
            featureSet = tw_extractFeatureSetFromDB_local(db,frames,'pos',[3,4,8,9,14,17,19,21,26,28]);
        %----------------------------------------------------------------------
        case 'e39' 
            featureSet = tw_extractFeatureSetFromDB_local(db,frames,'pos',[3,4,8,9,12,14,17,18,19,21,25,26,28]);
        %----------------------------------------------------------------------
        case 'e15_45'
            data        = tw_extractFeatureSetFromDB(db,'e15');
            offset      = [-10 -5 0 5 10];
            featureSet  = stackData(data,offset,db.motStartIDs);
            featureSet  = featureSet(:,frames);
        %----------------------------------------------------------------------
        case 'e15_flex'
            data        = tw_extractFeatureSetFromDB(db,'e15');
            offset      = frame_offsets;
            featureSet  = stackData(data,offset,db.motStartIDs);
            featureSet  = featureSet(:,frames);
        %----------------------------------------------------------------------
        case 'e15_75'
            data        = tw_extractFeatureSetFromDB(db,'e15');
            offset      = [-5 -3 0 3 5];
            featureSet  = stackData(data,offset,db.motStartIDs);
            featureSet  = featureSet(:,frames);
        %----------------------------------------------------------------------
        case 'e30_90'
            data        = tw_extractFeatureSetFromDB(db,'e30');
            offset      = [-5 0 5];
            featureSet  = stackData(data,offset,db.motStartIDs);
            featureSet  = featureSet(:,frames);
        %----------------------------------------------------------------------
        case 'e30_flex'
            data        = tw_extractFeatureSetFromDB(db,'e30');
            offset      = frame_offsets;
            featureSet  = stackData(data,offset,db.motStartIDs);
            featureSet  = featureSet(:,frames);
        %----------------------------------------------------------------------
        case 'e27_81'
            data        = tw_extractFeatureSetFromDB(db,'e27');
            offset      = [-5 0 5];
            featureSet  = stackData(data,offset,db.motStartIDs);
            featureSet  = featureSet(:,frames);
    end
end

function [featureSet,db] = tw_extractFeatureSetFromDB_local(db, frames, F, varargin)
        
    %----------------------------------------------------------------------
    % Generalised feature sets
    %----------------------------------------------------------------------
    
    switch lower(F)
        case {'pos','vel','acc','euler','quat','posorig','velorig'}
            switch nargin
                case 4
                    % The jointIDs identify the joints which are used in
                    % the feature set
                    jointIDs    = varargin{1};
                otherwise
                    error('Wrong number of argins');
            end
            
            matrixIDs = jointIDsToMatrixIndices(jointIDs);

            if isfield(db,F)
                switch F
                    case {'pos','vel','acc','posOrig','velOrig'}
                        featureSet = cast(db.(F)(matrixIDs.pos, frames),'double');
                end

            else
                error('Field not existent in db!');
            end
        %----------------------------------------------------------------------
        otherwise
            error('Unknown feature set!');
    end
end
