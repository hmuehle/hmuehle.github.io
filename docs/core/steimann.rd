<?xml version="1.0" encoding="UTF-8"?>
<roles:RoleModel xmi:version="2.0" xmlns:xmi="http://www.omg.org/XMI" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:roles="http://net.core.editor/1.0">
  <baseTypes name="SuperBase">
    <attributes name="bAtt1" type="String"/>
  </baseTypes>
  <baseTypes name="SubBase">
    <attributes name="bAtt2" type="String"/>
  </baseTypes>
  <edges xsi:type="roles:BaseInheritanceEdge" source="//@baseTypes.1" target="//@baseTypes.0"/>
  <edges xsi:type="roles:RoleInheritanceEdge" source="//@collaborations.0/@roleTypes.1" target="//@collaborations.0/@roleTypes.0"/>
  <edges xsi:type="roles:RolePlayEdge" source="//@collaborations.0/@roleTypes.1" target="//@baseTypes.0"/>
  <collaborations name="Roles">
    <roleTypes name="SuperRole">
      <attributes name="rAtt1" type="String"/>
    </roleTypes>
    <roleTypes name="SubRole">
      <attributes name="rAtt2" type="String"/>
    </roleTypes>
  </collaborations>
</roles:RoleModel>
