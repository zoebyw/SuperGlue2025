import React, { useState, useEffect, useRef } from 'react';
import { Button, message, Tabs } from 'antd';
import { useNavigate } from 'react-router-dom';
import '../styles/main.css';
import turnBackIcon from '../assets/back-arrow.png';

const MoleculeIndex = () => {
  const [ketcher, setKetcher] = useState(null);
  const [currentSmiles, setCurrentSmiles] = useState('');
  const [ketcherSmiles, setKetcherSmiles] = useState('');
  const [ketcherMolfile, setKetcherMolfile] = useState('');
  const iframeRef = useRef(null);
  const navigate = useNavigate();

  useEffect(() => {
    const initializeKetcher = () => {
      const ketcherFrame = document.getElementById('idKetcher');
      if (!ketcherFrame) {
        message.error('Failed to locate the iframe element.');
        return;
      }

      const ketcherInstance = ketcherFrame.contentWindow?.ketcher || document.frames['idKetcher']?.window?.ketcher;

      if (ketcherInstance) {
        setKetcher(ketcherInstance);
        message.success('Ketcher initialized successfully.');
      } else {
        message.error('Failed to initialize Ketcher.');
      }
    };

    const timer = setTimeout(() => {
      initializeKetcher();
    }, 500);

    return () => clearTimeout(timer);
  }, []);

  const applySmiles = () => {
    if (!ketcher) {
      message.error('Ketcher instance not initialized.');
      return;
    }
    try {
      ketcher.setMolecule(currentSmiles);
      message.success('SMILES applied successfully.');
    } catch (error) {
      console.error('Error applying SMILES to Ketcher:', error);
      message.error('Failed to apply SMILES.');
    }
  };

  const getSmiles = async () => {
    if (!ketcher) {
      message.error('Ketcher instance not initialized.');
      return;
    }
    try {
      const smiles = await ketcher.getSmiles();
      setKetcherSmiles(smiles);
      message.success(`SMILES: ${smiles}`);
    } catch (error) {
      console.error('Error fetching SMILES from Ketcher:', error);
      message.error('Failed to get SMILES.');
    }
  };

  const getMolfile = async () => {
    if (!ketcher) {
      message.error('Ketcher instance not initialized.');
      return;
    }
    try {
      const molfile = await ketcher.getMolfile();
      setKetcherMolfile(molfile);
      message.success('Molfile fetched successfully.');
    } catch (error) {
      console.error('Error fetching Molfile from Ketcher:', error);
      message.error('Failed to get Molfile.');
    }
  };

  // return button handler
  const handleBack = () => {
    const confirmExit = window.confirm('Do you want to save your changes before exiting?');
    if (confirmExit) {
      // save changes
      console.log('Changes saved');
    }
    navigate(-1); // return to previous page
  };

  return (
    <div className="editor-container">
      <div className="editor-header">
        <img
          src={turnBackIcon}
          alt="Return"
          className="back-icon"
          onClick={handleBack}
        />
        <h2>Chemical Canvas</h2>
      </div>
      <div className="editor-main">
        <div className="editor-left">
          <Tabs defaultActiveKey="ketcher" type="card">
            <Tabs.TabPane tab="Ketcher" key="ketcher">
              <iframe
                id="idKetcher"
                ref={iframeRef}
                src="./standalone/index.html"
                className="ketcher-frame"
              />
            </Tabs.TabPane>
          </Tabs>
        </div>
        <div className="editor-right">
          <div className="input-group">
            <input
              className="input-search"
              placeholder="Input SMILES text"
              onChange={(e) => setCurrentSmiles(e.target.value)}
              value={currentSmiles}
            />
            <Button onClick={applySmiles} type="primary">
              Apply
            </Button>
          </div>
          <div className="output-group">
            <Button onClick={getSmiles} type="primary">
              Get SMILES
            </Button>
            <Button onClick={getMolfile} type="primary">
              Get Molfile
            </Button>
            <p>SMILES: {ketcherSmiles}</p>
            <p>Molfile: {ketcherMolfile}</p>
          </div>
        </div>
      </div>
    </div>
  );
};

export default MoleculeIndex;
