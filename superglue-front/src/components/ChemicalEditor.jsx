import React, { useState, useEffect, useRef } from 'react';
import { Button, message, Tabs, Layout, Menu, Card, List, Typography, Table, Space, Divider } from 'antd';
import { useNavigate, useLocation, useParams } from 'react-router-dom';
import {
  EditOutlined,
  FileSearchOutlined,
  CalculatorOutlined,
  ExportOutlined,
  CommentOutlined,
  ArrowLeftOutlined,
  SearchOutlined,
  TableOutlined
} from '@ant-design/icons';
import SimilaritySearch from './SimilaritySearch';
import '../styles/main.css';

const { Header, Sider, Content } = Layout;
const { Title, Text } = Typography;
const { TabPane } = Tabs;

const MoleculeIndex = () => {
  // Get the URL parameter
  const { id } = useParams();
  const location = useLocation();
  const navigate = useNavigate();
  const iframeRef = useRef(null);

  // Use the URL id parameter
  const moleculeIdFromParams = id ? parseInt(id, 10) : null;

  const [ketcher, setKetcher] = useState(null);
  const [currentSmiles, setCurrentSmiles] = useState('');
  const [ketcherSmiles, setKetcherSmiles] = useState('');
  const [ketcherMolfile, setKetcherMolfile] = useState('');
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const [selectedTab, setSelectedTab] = useState('annotate');

 // state for similarity search
  const [similaritySearchVisible, setSimilaritySearchVisible] = useState(false);
  const [showResultsTable, setShowResultsTable] = useState(false);
  // const [searchResults, setSearchResults] = useState([]);
  const [similarityResults, setSimilarityResults] = useState([]);
  const [isSearching, setIsSearching] = useState(false);
  const [searchError, setSearchError] = useState('');
  const [similarityMethod, setSimilarityMethod] = useState('tanimoto');
  const [searchQuery, setSearchQuery] = useState('');
  const [searchMethod, setSearchMethod] = useState('');


  // Get dynamic molecule property from CSV
  const [moleculeProperties, setMoleculeProperties] = useState({});
  const [propertyKeys, setPropertyKeys] = useState([]);

  // Get data from last webpage
  const smilesFromCSV = location.state?.smiles;
  const moleculeId = location.state?.moleculeId || moleculeIdFromParams;
  const moleculeName = location.state?.moleculeName || `Compound-${moleculeIdFromParams}`;
  const fromCsv = location.state?.fromCsv;
  const sourceData = location.state?.sourceData; // properties
  const idColumn = location.state?.idColumn || "cmpd_id";
  const smilesColumn = location.state?.smilesColumn || "SMILES";
  const filename = location.state?.filename;

  // Logging for debugging
  useEffect(() => {
    console.log("URL Parameter ID:", id);
    console.log("Location state:", location.state);
    console.log("SMILES from CSV:", smilesFromCSV);
    console.log("ID column:", idColumn);
    console.log("SMILES column:", smilesColumn);
    console.log("Source Data:", sourceData);
  }, [location, id]);

  // Initialize smiles
  useEffect(() => {
    if (smilesFromCSV) {
      setCurrentSmiles(smilesFromCSV);
    }
  }, [smilesFromCSV]);

  // Show property
  useEffect(() => {
    if (sourceData) {
      // exclude property
      const basicFields = [idColumn, smilesColumn, 'id', 'smiles', 'structure'];

      // extract property from sourceData, exclude id and smiles
      const properties = {};
      const propertyNames = [];

      Object.entries(sourceData).forEach(([key, value]) => {
        // check if it is id or smiles
        const isBasicField = basicFields.some(field =>
          key.toLowerCase() === field.toLowerCase() ||
          (field !== idColumn && field !== smilesColumn && (
            key.toLowerCase().includes('id') ||
            key.toLowerCase().includes('name') ||
            key.toLowerCase().includes('smiles')
          ))
        );

        if (!isBasicField && value !== undefined && value !== null && value !== '') {
          properties[key] = value;
          propertyNames.push(key);
        }
      });

      setMoleculeProperties(properties);
      setPropertyKeys(propertyNames);

      console.log("Extracted properties:", properties);
    }
  }, [sourceData, idColumn, smilesColumn]);

  // Initialize Ketcher
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
        if (smilesFromCSV) {
          setTimeout(() => {
            try {
              ketcherInstance.setMolecule(smilesFromCSV);
              message.success('Molecule loaded successfully');
            } catch (error) {
              console.error('Error applying SMILES from CSV:', error);
              message.error('Failed to load molecule structure.');
            }
          }, 800); // give some time for ketcher to visualize
        }
      } else {
        message.error('Failed to initialize Ketcher.');
      }
    };

    const timer = setTimeout(() => {
      initializeKetcher();
    }, 500);

    return () => clearTimeout(timer);
  }, []);

  // Visualize smiles in Ketcher
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

  // Get smiles of molecule
  const getSmiles = async () => {
    if (!ketcher) {
      message.error('Ketcher instance not initialized.');
      return;
    }
    try {
      const smiles = await ketcher.getSmiles();
      setKetcherSmiles(smiles);
      setCurrentSmiles(smiles); // Update current SMILES with the one from Ketcher
      message.success(`SMILES: ${smiles}`);
    } catch (error) {
      console.error('Error fetching SMILES from Ketcher:', error);
      message.error('Failed to get SMILES.');
    }
  };

  // Get molfile of molecule
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

  // Return to last page
  const handleBack = () => {
    const confirmExit = window.confirm('Do you want to save your changes before exiting?');
    if (confirmExit) {
      // Save changes logic would go here
      console.log('Changes saved');
    }
    navigate(-1);
  };

  // Handle sidebar actions
  const handleSidebarAction = (action) => {
    setSelectedTab(action);
    message.info(`${action.charAt(0).toUpperCase() + action.slice(1)} action selected`);

    // Handle special functionality
    if (action === 'compute') {
      getSmiles();
    }

    // Handle similarity search panel
    if (action === 'similarity') {
      // First update current SMILES from Ketcher
      getSmiles().then(() => {
        setSimilaritySearchVisible(true);
      });
    } else {
      setSimilaritySearchVisible(false);
      setShowResultsTable(false);
    }
  };

  // Handle similarity search results
  const handleSimilarityResults = (results, method) => {
    setSimilarityResults(results);
    setShowResultsTable(true);
    setSearchMethod(method);

    getSmiles().then(smiles => {
    setSearchQuery(smiles || currentSmiles);
  });
  };

  // Format property values
  const formatPropertyValue = (value) => {
    if (typeof value === 'number') {
      if (Number.isInteger(value)) {
        return value.toString();
      } else {
        return value.toFixed(2);
      }
    } else if (typeof value === 'string') {
      const numValue = parseFloat(value);
      if (!isNaN(numValue) && value.trim() !== '') {
        return numValue.toFixed(2);
      }
      return value;
    }
    return String(value);
  };

  // Get property card style based on value
  const getPropertyCardStyle = (key, value) => {
    const numValue = parseFloat(value);
    if (isNaN(numValue)) {
      return {};
    }

    if (key.toLowerCase().includes('logp') || key.toLowerCase().includes('alogp')) {
      if (numValue > 5) return { backgroundColor: '#ffecec' };
      if (numValue < 0) return { backgroundColor: '#ecffec' };
    } else if (key.toLowerCase().includes('weight') || key.toLowerCase().includes('mw')) {
      if (numValue > 500) return { backgroundColor: '#f0f0ff' };
      if (numValue < 200) return { backgroundColor: '#fffff0' };
    } else if (key.toLowerCase().includes('tpsa') || key.toLowerCase().includes('psa')) {
      if (numValue > 140) return { backgroundColor: '#e6f7ff' };
      if (numValue < 60) return { backgroundColor: '#fffbe6' };
    }

    return {};
  };

  // Results table columns
  const resultColumns = [
    {
    title: 'Cmpd Id',
    dataIndex: 'cmpd_id',
    key: 'cmpd_id',
  },
  {
    title: 'Similarity',
    dataIndex: 'similarity',
    key: 'similarity',
    sorter: (a, b) => a.similarity - b.similarity,
    render: value => (value * 100).toFixed(1) + '%',
    defaultSortOrder: 'descend',
  },

  {
    title: 'Binary Occ',
    dataIndex: 'binary_occ',
    key: 'binary_occ',
    sorter: (a, b) => a.binary_occ - b.binary_occ,
    render: value => value?.toFixed(2) || '-',
  },
  {
    title: 'Cont Occ',
    dataIndex: 'cont_occ',
    key: 'cont_occ',
    sorter: (a, b) => a.cont_occ - b.cont_occ,
    render: value => value?.toFixed(2) || '-',
  },
  {
    title: 'Low Gsh Prob',
    dataIndex: 'low_gsh_prob',
    key: 'low_gsh_prob',
    sorter: (a, b) => a.low_gsh_prob - b.low_gsh_prob,
    render: value => value?.toFixed(2) || '-',
  },
  {
    title: 'Med Gsh Prob',
    dataIndex: 'med_gsh_prob',
    key: 'med_gsh_prob',
    sorter: (a, b) => a.med_gsh_prob - b.med_gsh_prob,
    render: value => value?.toFixed(2) || '-',
  },
  {
    title: 'High Gsh Prob',
    dataIndex: 'high_gsh_prob',
    key: 'high_gsh_prob',
    sorter: (a, b) => a.high_gsh_prob - b.high_gsh_prob,
    render: value => value?.toFixed(2) || '-',
  },
  {
    title: 'Selectivity',
    dataIndex: 'selectivity',
    key: 'selectivity',
    sorter: (a, b) => a.selectivity - b.selectivity,
    render: value => value?.toFixed(2) || '-',
  },
  {
    title: 'Actions',
    key: 'actions',
    render: (_, record) => (
      <Space>
        <Button size="small" onClick={() => {
          if (ketcher) {
            try {
              ketcher.setMolecule(record.smiles);
            } catch (error) {
              console.error('Error loading molecule:', error);
            }
          }
        }}>
          View
        </Button>
      </Space>
    ),
  },
];

  const ketcherPath = window.location.origin + '/standalone/index.html';

  return (
    <Layout style={{ minHeight: '100vh' }}>
      {/* Left Sidebar */}
      <Sider width={200} theme="light" style={{ boxShadow: '2px 0 8px rgba(0,0,0,0.15)' }}>
        <div className="logo" style={{ height: '64px', display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
          <Title level={4} style={{ margin: 0, color: '#1890ff' }}>Chemical Canvas</Title>
        </div>

        <div style={{ display: 'flex', flexDirection: 'column', height: 'calc(100% - 64px)' }}>
          <Menu
            mode="vertical"
            selectedKeys={[selectedTab]}
            style={{ borderRight: 0 }}
          >
            <Menu.Item key="modify" icon={<EditOutlined />} onClick={() => handleSidebarAction('modify')}>
              Modify
            </Menu.Item>
            <Menu.Item key="annotate" icon={<CommentOutlined />} onClick={() => handleSidebarAction('annotate')}>
              Annotate
            </Menu.Item>
            <Menu.Item key="similarity" icon={<FileSearchOutlined />} onClick={() => handleSidebarAction('similarity')}>
              Similarity Search
            </Menu.Item>
            <Menu.Item key="compute" icon={<CalculatorOutlined />} onClick={() => handleSidebarAction('compute')}>
              Compute
            </Menu.Item>
            <Menu.Item key="export" icon={<ExportOutlined />} onClick={() => handleSidebarAction('export')}>
              Export
            </Menu.Item>
          </Menu>

          <div style={{ marginTop: 'auto', padding: '20px', textAlign: 'center' }}>
            <Button icon={<ArrowLeftOutlined />} onClick={handleBack}>
              Back
            </Button>
          </div>
        </div>
      </Sider>

      <Layout>
        <Header style={{ background: '#fff', padding: '0 16px', display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <div>
            <Text strong>Current Mode: </Text>
            <Text>{selectedTab.charAt(0).toUpperCase() + selectedTab.slice(1)}</Text>
          </div>
          {showResultsTable && selectedTab === 'similarity' && (
            <div>
              <Button
                type="primary"
                icon={<TableOutlined />}
                onClick={() => setShowResultsTable(!showResultsTable)}
              >
                {showResultsTable ? 'Hide Results' : 'Show Results'}
              </Button>
            </div>
          )}
          <div>
            <Text strong>Molecule: </Text>
            <Text>{moleculeName}</Text>
          </div>
        </Header>

        <Layout style={{ padding: '24px' }}>
          {/* Main content layout with conditional height based on results visibility */}
          <Layout>
            {/* Ketcher editor area */}
            <Content
              style={{
                background: '#fff',
                padding: '12px',
                width: '100%',
                height: showResultsTable ? 'calc(60vh - 112px)' : 'calc(100vh - 112px)',
                marginBottom: showResultsTable ? '16px' : '0',
                transition: 'height 0.3s ease'
              }}
            >
              <iframe
                id="idKetcher"
                ref={iframeRef}
                src={ketcherPath}
                style={{ width: '100%', height: '100%', border: 'none' }}
              />
            </Content>

            {/* Similarity Results Table (conditionally rendered) */}
            {showResultsTable && (
              <Content style={{ background: '#fff', padding: '12px', width: '100%', marginBottom: '16px' }}>
                <div style={{ marginBottom: '12px' }}>
                  <div style={{ marginBottom: '8px' }}>
                    <Text strong>Query Smiles: </Text>
                    <Text>{searchQuery}</Text>
                  </div>
                  <div style={{ marginBottom: '8px' }}>
                    <Text strong>Similarity Metric: </Text>
                    <Text>{searchMethod.charAt(0).toUpperCase() + searchMethod.slice(1)}</Text>
                  </div>
                  <Text strong>Search Results</Text>
                </div>

                <Table
                  columns={resultColumns}
                  dataSource={similarityResults}
                  size="small"
                  pagination={{
                    pageSize: 5,
                    showSizeChanger: true,
                    pageSizeOptions: ['5', '10', '20']
                  }}
                  scroll={{ y: 240 }}
                />
              </Content>
            )}

          </Layout>

          {/* Right sidebar with properties */}
          <Sider width="30%" style={{ background: '#f0f2f5', marginLeft: '16px' }}>
            {/* Basic Information card */}
            <Card title="Basic Information" style={{ marginBottom: '16px' }}>
              <List size="small">
                <List.Item>
                  <Text strong style={{ width: '40%' }}>ID:</Text>
                  <Text>{moleculeName}</Text>
                </List.Item>
                {sourceData && sourceData[smilesColumn] && (
                  <List.Item>
                    <Text strong style={{ width: '40%' }}>SMILES:</Text>
                    <Text style={{ wordBreak: 'break-all' }}>
                      {sourceData[smilesColumn].length > 50
                        ? sourceData[smilesColumn].substring(0, 50) + '...'
                        : sourceData[smilesColumn]
                      }
                    </Text>
                  </List.Item>
                )}
                {filename && (
                  <List.Item>
                    <Text strong style={{ width: '40%' }}>Source:</Text>
                    <Text>{filename}</Text>
                  </List.Item>
                )}
              </List>
            </Card>

            {/* Properties from CSV */}
            {propertyKeys.length > 0 ? (
              <Card title="Properties from CSV" style={{ marginBottom: '16px' }}>
                <List
                  size="small"
                  itemLayout="horizontal"
                  dataSource={propertyKeys}
                  renderItem={key => {
                    const value = moleculeProperties[key];
                    return (
                      <List.Item style={getPropertyCardStyle(key, value)}>
                        <Text strong style={{ width: '50%' }}>{key}:</Text>
                        <Text>{formatPropertyValue(value)}</Text>
                      </List.Item>
                    );
                  }}
                />
              </Card>
            ) : (
              fromCsv && (
                <Card title="Properties from CSV" style={{ marginBottom: '16px' }}>
                  <Text type="secondary">No additional properties found in CSV.</Text>
                </Card>
              )
            )}

            {/* Mode-specific content cards */}
            {selectedTab === 'annotate' && (
              <Card title="Annotations">
                <textarea
                  placeholder="Add your annotations here..."
                  style={{ width: '100%', height: '100px', padding: '8px', borderRadius: '4px', border: '1px solid #d9d9d9' }}
                />
                <Button type="primary" size="small" style={{ marginTop: '8px' }}>
                  Save Annotations
                </Button>
              </Card>
            )}

            {selectedTab === 'compute' && (
              <Card title="Calculated Properties">
                <div style={{ marginBottom: '10px' }}>
                  <Button type="primary" size="small" onClick={getSmiles}>
                    Calculate Properties
                  </Button>
                </div>
                {ketcherSmiles && (
                  <div>
                    <Text strong>SMILES:</Text>
                    <div style={{ wordBreak: 'break-all', margin: '5px 0 10px' }}>
                      <Text>{ketcherSmiles}</Text>
                    </div>
                  </div>
                )}
              </Card>
            )}

            {selectedTab === 'export' && (
              <Card title="Export Options">
                <Button style={{ margin: '5px' }} onClick={getSmiles}>Get SMILES</Button>
                <Button style={{ margin: '5px' }} onClick={getMolfile}>Get Molfile</Button>
                {ketcherSmiles && (
                  <div style={{ marginTop: '10px' }}>
                    <Text strong>SMILES:</Text>
                    <div style={{
                      padding: '5px',
                      background: '#f5f5f5',
                      borderRadius: '4px',
                      marginTop: '5px',
                      wordBreak: 'break-all'
                    }}>
                      {ketcherSmiles}
                    </div>
                  </div>
                )}
                {ketcherMolfile && (
                  <div style={{ marginTop: '10px' }}>
                    <Text strong>Molfile:</Text>
                    <div style={{
                      padding: '5px',
                      background: '#f5f5f5',
                      borderRadius: '4px',
                      marginTop: '5px',
                      height: '100px',
                      overflow: 'auto',
                      whiteSpace: 'pre-wrap',
                      fontSize: '12px'
                    }}>
                      {ketcherMolfile}
                    </div>
                  </div>
                )}
              </Card>
            )}
          </Sider>
        </Layout>
      </Layout>

      {/* Similarity Search Modal */}
      <SimilaritySearch
        currentSmiles={currentSmiles}
        currentId={moleculeName || moleculeId || `Compound-${moleculeIdFromParams}`}
        filename={filename}
        visible={similaritySearchVisible}
        onClose={() => setSimilaritySearchVisible(false)}
        onResultsFound={handleSimilarityResults}
        ketcher={ketcher}
        moleculeProperties={moleculeProperties}
      />
    </Layout>
  );
};

export default MoleculeIndex;